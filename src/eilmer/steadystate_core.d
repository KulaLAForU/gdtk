/** steadystate_core.d
 * Core set of functions used in the Newton-Krylov updates for steady-state convergence.
 *
 * Author: Rowan G.
 * Date: 2016-10-09
 */

module steadystate_core;

import core.stdc.stdlib : exit;
import std.stdio;
import std.file;
import std.format;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.string;
import std.array;
import std.math;
import std.datetime;

import nm.smla;
import nm.bbla;
import nm.complex;
import nm.number;

import special_block_init;
import fluidblock;
import sfluidblock;
import globaldata;
import globalconfig;
import simcore;
import fvcore;
import fileutil;
import user_defined_source_terms;
import conservedquantities;
import postprocess : readTimesFile;
import loads;
import shape_sensitivity_core : sss_preconditioner_initialisation, sss_preconditioner;

static int fnCount = 0;
static shared bool with_k_omega;

// Module-local, global memory arrays and matrices
number[] g0;
number[] g1;
number[] h;
number[] hR;
Matrix!number H0;
Matrix!number H1;
Matrix!number Gamma;
Matrix!number Q0;
Matrix!number Q1;

struct RestartInfo {
    double pseudoSimTime;
    double dt;
    int step;
    double globalResidual;
    ConservedQuantities residuals;

    this(int n_species, int n_modes)
    {
        residuals = new ConservedQuantities(n_species, n_modes);
    }
}

void extractRestartInfoFromTimesFile(string jobName, ref RestartInfo[] times)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;
    size_t OMEGA = omegaIdx;

    auto gmodel = GlobalConfig.gmodel_master;
    RestartInfo restartInfo = RestartInfo(gmodel.n_species, gmodel.n_modes);
    // Start reading the times file, looking for the snapshot index
    auto timesFile = File(jobName ~ ".times");
    auto line = timesFile.readln().strip();
    while (line.length > 0) {
        if (line[0] != '#') {
            // Process a non-comment line.
            auto tokens = line.split();
            auto idx = to!int(tokens[0]);
            restartInfo.pseudoSimTime = to!double(tokens[1]);
            restartInfo.dt = to!double(tokens[2]);
            restartInfo.step = to!int(tokens[3]);
            restartInfo.globalResidual = to!double(tokens[4]);
            restartInfo.residuals.mass = to!double(tokens[5+MASS]);
            restartInfo.residuals.momentum.refx = to!double(tokens[5+X_MOM]);
            restartInfo.residuals.momentum.refy = to!double(tokens[5+Y_MOM]);
            if ( GlobalConfig.dimensions == 3 ) 
                restartInfo.residuals.momentum.refz = to!double(tokens[5+Z_MOM]);
            restartInfo.residuals.total_energy = to!double(tokens[5+TOT_ENERGY]);
            if ( GlobalConfig.turbulence_model == TurbulenceModel.k_omega ) {
                restartInfo.residuals.tke = to!double(tokens[5+TKE]);
                restartInfo.residuals.omega = to!double(tokens[5+OMEGA]);
            }
            times ~= restartInfo;
        }
        line = timesFile.readln().strip();
    }
    timesFile.close();
    return;
} 


void iterate_to_steady_state(int snapshotStart, int maxCPUs)
{
    auto wallClockStart = Clock.currTime();
    bool withPTC = true;
    string jobName = GlobalConfig.base_file_name;
    int nsteps = GlobalConfig.sssOptions.nTotalSteps;
    int nRestarts;
    int maxNumberAttempts = GlobalConfig.sssOptions.maxNumberAttempts;
    double relGlobalResidReduction = GlobalConfig.sssOptions.stopOnRelGlobalResid;
    double absGlobalResidReduction = GlobalConfig.sssOptions.stopOnAbsGlobalResid;

    // Settings for start-up phase
    double cfl0 = GlobalConfig.sssOptions.cfl0;
    double tau0 = GlobalConfig.sssOptions.tau0;
    double eta0 = GlobalConfig.sssOptions.eta0;
    double sigma0 = GlobalConfig.sssOptions.sigma0;
    // Settings for inexact Newton phase
    double cfl1 = GlobalConfig.sssOptions.cfl1;
    double tau1 = GlobalConfig.sssOptions.tau1;
    double sigma1 = GlobalConfig.sssOptions.sigma1;
    EtaStrategy etaStrategy = GlobalConfig.sssOptions.etaStrategy;
    double eta1 = GlobalConfig.sssOptions.eta1;
    double eta1_max = GlobalConfig.sssOptions.eta1_max;
    double eta1_min = GlobalConfig.sssOptions.eta1_min;
    double etaRatioPerStep = GlobalConfig.sssOptions.etaRatioPerStep;
    double gamma = GlobalConfig.sssOptions.gamma;
    double alpha = GlobalConfig.sssOptions.alpha;

    int interpOrderSave = GlobalConfig.interpolation_order;
    
    int n_species = GlobalConfig.gmodel_master.n_species();
    int n_modes = GlobalConfig.gmodel_master.n_modes();
    ConservedQuantities maxResiduals = new ConservedQuantities(n_species, n_modes);
    ConservedQuantities currResiduals = new ConservedQuantities(n_species, n_modes);

    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;
    size_t OMEGA = omegaIdx;

    shared double dt;
    double dtTrial, etaTrial;
    shared bool failedAttempt;
    double pseudoSimTime = 0.0;
    double normOld, normNew;
    int snapshotsCount = GlobalConfig.sssOptions.snapshotsCount;
    int nTotalSnapshots = GlobalConfig.sssOptions.nTotalSnapshots;
    int nWrittenSnapshots = 0;
    int writeDiagnosticsCount = GlobalConfig.sssOptions.writeDiagnosticsCount;
    int writeLoadsCount = GlobalConfig.sssOptions.writeLoadsCount;

    int startStep;
    int nPreSteps = GlobalConfig.sssOptions.nPreSteps;
    int nStartUpSteps = GlobalConfig.sssOptions.nStartUpSteps;
    bool inexactNewtonPhase = false;

    // No need to have more task threads than blocks
    auto nThreadsInPool = min(maxCPUs-1, GlobalConfig.nFluidBlocks-1);
    defaultPoolThreads(nThreadsInPool);
    writeln("Running with ", nThreadsInPool+1, " threads."); // +1 for main thread.

    double normRef = 0.0;
    bool residualsUpToDate = false;
    bool finalStep = false;
    bool usePreconditioner = GlobalConfig.sssOptions.usePreconditioner;
    if (usePreconditioner) {
        writeln("Initialising memory.");
        foreach (blk; parallel(localFluidBlocks,1)) {
            sss_preconditioner_initialisation(blk, nConserved); 
        }
    }
    // Set usePreconditioner to false for pre-steps AND first-order steps.
    usePreconditioner = false;

    // We only do a pre-step phase if we are starting from scratch.
    if ( snapshotStart == 0 ) {
        dt = determine_initial_dt(cfl0);
        // The initial residual is usually a poor choice for basing decisions about how the
        // residual is dropping, particularly when a constant initial condition is given.
        // A constant initial condition gives a zero residual everywhere in the interior
        // of the domain. Therefore, the only residual can come from boundary condition influences.
        // What we do here is use some early steps (I've called them "pre"-steps) to
        // allow the boundary conditions to exert some influence into the interior.
        // We'll find the max residual in that start-up process and use that as our initial value
        // for the "max" residual. Our relative residuals will then be based on that.
        // We can use a fixed timestep (typically small) and low-order reconstruction in this
        // pre-step phase.

        writeln("Begin pre-iterations to establish sensible max residual.");
        foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(1);
        foreach ( preStep; -nPreSteps .. 0 ) {
            foreach (attempt; 0 .. maxNumberAttempts) {
                try {
                    rpcGMRES_solve(pseudoSimTime, dt, eta0, sigma0, usePreconditioner, normOld, nRestarts);
                }
                catch (FlowSolverException e) {
                    writefln("Failed when attempting GMRES solve in pre-steps.");
                    writefln("attempt %d: dt= %e", attempt, dt);
                    failedAttempt = true;
                    dt = 0.1*dt;
                    break;
                }
                foreach (blk; parallel(localFluidBlocks,1)) {
                    bool local_with_k_omega = with_k_omega;
                    int cellCount = 0;
                    foreach (cell; blk.cells) {
                        cell.U[1].copy_values_from(cell.U[0]);
                        cell.U[1].mass = cell.U[0].mass + blk.dU[cellCount+MASS];
                        cell.U[1].momentum.refx = cell.U[0].momentum.x + blk.dU[cellCount+X_MOM];
                        cell.U[1].momentum.refy = cell.U[0].momentum.y + blk.dU[cellCount+Y_MOM];
                        if ( blk.myConfig.dimensions == 3 ) 
                            cell.U[1].momentum.refz = cell.U[0].momentum.z + blk.dU[cellCount+Z_MOM];
                        cell.U[1].total_energy = cell.U[0].total_energy + blk.dU[cellCount+TOT_ENERGY];
                        if ( local_with_k_omega ) {
                            cell.U[1].tke = cell.U[0].tke + blk.dU[cellCount+TKE];
                            cell.U[1].omega = cell.U[0].omega + blk.dU[cellCount+OMEGA];
                        }
                        try {
                            cell.decode_conserved(0, 1, 0.0);
                        }
                        catch (FlowSolverException e) {
                            writefln("Failed to provide sensible update.");
                            writefln("attempt %d: dt= %e", attempt, dt);
                            failedAttempt = true;
                            dt = 0.1*dt;
                            break;
                        }
                        cellCount += nConserved;
                    }
                }

                if ( failedAttempt )
                    continue;

                // If we get here, things are good. Put flow state into U[0]
                // ready for next iteration.
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (cell; blk.cells) {
                        swap(cell.U[0], cell.U[1]);
                    }
                }
                // If we got here, we can break the attempts loop.
                break;
            }
            if ( failedAttempt ) {
                writefln("Pre-step failed: %d", step);
                writeln("Bailing out!");
                exit(1);
            }
        
            writefln("PRE-STEP  %d  ::  dt= %.3e  global-norm= %.12e", preStep, dt, normOld); 
            if ( normOld > normRef ) {
                normRef = normOld;
                max_residuals(maxResiduals);
            }
        }

        writeln("Pre-step phase complete.");
    
        if ( nPreSteps <= 0 ) {
            // Take initial residual as max residual
            evalRHS(0.0, 0);
            max_residuals(maxResiduals);
            foreach (blk; parallel(localFluidBlocks, 1)) {
                bool local_with_k_omega = with_k_omega;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    blk.FU[cellCount+MASS] = -cell.dUdt[0].mass;
                    blk.FU[cellCount+X_MOM] = -cell.dUdt[0].momentum.x;
                    blk.FU[cellCount+Y_MOM] = -cell.dUdt[0].momentum.y;
                    if ( GlobalConfig.dimensions == 3 )
                        blk.FU[cellCount+Z_MOM] = -cell.dUdt[0].momentum.z;
                    blk.FU[cellCount+TOT_ENERGY] = -cell.dUdt[0].total_energy;
                    if ( local_with_k_omega ) {
                        blk.FU[cellCount+TKE] = -cell.dUdt[0].tke;
                        blk.FU[cellCount+OMEGA] = -cell.dUdt[0].omega;
                    }
                    cellCount += nConserved;
                }
            }
            mixin(norm2_over_blocks("normRef", "FU"));
        }
        writeln("Reference residuals are established as:");
        writefln("GLOBAL:         %.12e", normRef);
        writefln("MASS:           %.12e", maxResiduals.mass);
        writefln("X-MOMENTUM:     %.12e", maxResiduals.momentum.x);
        writefln("Y-MOMENTUM:     %.12e", maxResiduals.momentum.y);
        if ( GlobalConfig.dimensions == 3 )
            writefln("Z-MOMENTUM:     %.12e", maxResiduals.momentum.z);
        writefln("ENERGY:         %.12e", maxResiduals.total_energy);
        if ( with_k_omega ) {
            writefln("TKE:            %.12e", maxResiduals.tke);
            writefln("OMEGA:          %.12e", maxResiduals.omega);
        }
        
        string refResidFname = jobName ~ "-ref-residuals.saved";
        auto refResid = File(refResidFname, "w");
        if ( GlobalConfig.dimensions == 2 ) {
            refResid.writef("%.18e %.18e %.18e %.18e %.18e",
                            normRef, maxResiduals.mass, maxResiduals.momentum.x,
                            maxResiduals.momentum.y, maxResiduals.total_energy);
        }
        else {
            refResid.writef("%.18e %.18e %.18e %.18e %.18e %.18e",
                            normRef, maxResiduals.mass, maxResiduals.momentum.x,
                            maxResiduals.momentum.y, maxResiduals.momentum.z,
                            maxResiduals.total_energy);
        }
        if ( with_k_omega ) {
            refResid.writef(" %.18e %.18e", maxResiduals.tke, maxResiduals.omega);
        }
        refResid.write("\n");
        refResid.close();
    }

    RestartInfo[] times;

    if ( snapshotStart > 0 ) {
        extractRestartInfoFromTimesFile(jobName, times);
        normOld = times[snapshotStart].globalResidual;
        // We need to read in the reference residual values from a file.
        string refResidFname = jobName ~ "-ref-residuals.saved";
        auto refResid = File(refResidFname, "r");
        auto line = refResid.readln().strip();
        auto tokens = line.split();
        normRef = to!double(tokens[0]);
        maxResiduals.mass = to!double(tokens[1+MASS]);
        maxResiduals.momentum.refx = to!double(tokens[1+X_MOM]);
        maxResiduals.momentum.refy = to!double(tokens[1+Y_MOM]);
        if ( GlobalConfig.dimensions == 3 ) 
            maxResiduals.momentum.refz = to!double(tokens[1+Z_MOM]);
        maxResiduals.total_energy = to!double(tokens[1+TOT_ENERGY]);
        if ( with_k_omega ) {
            maxResiduals.tke = to!double(tokens[1+TKE]);
            maxResiduals.omega = to!double(tokens[1+OMEGA]);
        }
        
        // We also need to determine how many snapshots have already been written
        auto timesFile = File(jobName ~ ".times");
        line = timesFile.readln().strip();
        while (line.length > 0) {
            if (line[0] != '#') {
                nWrittenSnapshots++;
            }
            line = timesFile.readln().strip();
        }
        timesFile.close();
        nWrittenSnapshots--; // We don't count the initial solution as a written snapshot
    }

    double wallClockElapsed;
    double cfl;
    RestartInfo restartInfo;

    // We need to do some configuration based on whether we are starting from scratch,
    // or attempting to restart from an earlier snapshot.
    if ( snapshotStart == 0 ) {
        startStep = 1;
        restartInfo.pseudoSimTime = 0.0;
        restartInfo.dt = dt;
        restartInfo.step = 0;
        restartInfo.globalResidual = normRef;
        restartInfo.residuals = maxResiduals;
        times ~= restartInfo;
    }
    else {
        restartInfo = times[snapshotStart];
        dt = restartInfo.dt;
        startStep = restartInfo.step + 1;
        pseudoSimTime = restartInfo.pseudoSimTime;
        writefln("Restarting steps from step= %d", startStep);
        writefln("   pseudo-sim-time= %.6e dt= %.6e", pseudoSimTime, dt);  
    }


    auto residFname = "e4sss.diagnostics.dat";
    File fResid;
    if ( snapshotStart == 0 ) {
        // Open file for writing diagnostics
        fResid = File(residFname, "w");
        fResid.writeln("#  1: step");
        fResid.writeln("#  2: pseudo-time");
        fResid.writeln("#  3: dt");
        fResid.writeln("#  4: CFL");
        fResid.writeln("#  5: eta");
        fResid.writeln("#  6: nRestarts");
        fResid.writeln("#  7: nFnCalls");
        fResid.writeln("#  8: wall-clock, s");
        fResid.writeln("#  9: global-residual-abs");
        fResid.writeln("# 10: global-residual-rel");
        fResid.writefln("#  %02d: mass-abs", 11+2*MASS);
        fResid.writefln("# %02d: mass-rel", 11+2*MASS+1);
        fResid.writefln("# %02d: x-mom-abs", 11+2*X_MOM);
        fResid.writefln("# %02d: x-mom-rel", 11+2*X_MOM+1);
        fResid.writefln("# %02d: y-mom-abs", 11+2*Y_MOM);
        fResid.writefln("# %02d: y-mom-rel", 11+2*Y_MOM+1);
        if ( GlobalConfig.dimensions == 3 ) {
            fResid.writefln("# %02d: z-mom-abs", 11+2*Z_MOM);
            fResid.writefln("# %02d: z-mom-rel", 11+2*Z_MOM+1);
        }
        fResid.writefln("# %02d: energy-abs", 11+2*TOT_ENERGY);
        fResid.writefln("# %02d: energy-rel", 11+2*TOT_ENERGY+1);
        if ( with_k_omega ) {
            fResid.writefln("# %02d: tke-abs", 11+2*TKE);
            fResid.writefln("# %02d: tke-rel", 11+2*TKE+1);
            fResid.writefln("# %02d: omega-abs", 11+2*OMEGA);
            fResid.writefln("# %02d: omega-rel", 11+2*OMEGA+1);
        }

        fResid.close();
    }

    // Begin Newton steps
    double eta;
    double tau;
    double sigma;
    foreach (step; startStep .. nsteps+1) {
        if ( (step/GlobalConfig.control_count)*GlobalConfig.control_count == step ) {
            read_control_file(); // Reparse the time-step control parameters occasionally.
        }
        residualsUpToDate = false;
        if ( step <= nStartUpSteps ) {
            foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(1);
            eta = eta0;
            tau = tau0;
            sigma = sigma0;
            usePreconditioner = false;
        }
        else {
            foreach (blk; parallel(localFluidBlocks,1)) blk.set_interpolation_order(interpOrderSave);
            eta = eta1;
            tau = tau1;
            sigma = sigma1;
            usePreconditioner = GlobalConfig.sssOptions.usePreconditioner;
        }

        foreach (attempt; 0 .. maxNumberAttempts) {
            failedAttempt = false;
            try {
                rpcGMRES_solve(pseudoSimTime, dt, eta, sigma, usePreconditioner, normNew, nRestarts);
            }
            catch (FlowSolverException e) {
                writefln("Failed when attempting GMRES solve in main steps.");
                writefln("attempt %d: dt= %e", attempt, dt);
                failedAttempt = true;
                dt = 0.1*dt;
                continue;
            }
            foreach (blk; parallel(localFluidBlocks,1)) {
                bool local_with_k_omega = with_k_omega;
                int cellCount = 0;
                foreach (cell; blk.cells) {
                    cell.U[1].copy_values_from(cell.U[0]);
                    cell.U[1].mass = cell.U[0].mass + blk.dU[cellCount+MASS];
                    cell.U[1].momentum.refx = cell.U[0].momentum.x + blk.dU[cellCount+X_MOM];
                    cell.U[1].momentum.refy = cell.U[0].momentum.y + blk.dU[cellCount+Y_MOM];
                    if ( blk.myConfig.dimensions == 3 ) 
                        cell.U[1].momentum.refz = cell.U[0].momentum.z + blk.dU[cellCount+Z_MOM];
                    cell.U[1].total_energy = cell.U[0].total_energy + blk.dU[cellCount+TOT_ENERGY];
                    if ( local_with_k_omega ) {
                        cell.U[1].tke = cell.U[0].tke + blk.dU[cellCount+TKE];
                        cell.U[1].omega = cell.U[0].omega + blk.dU[cellCount+OMEGA];
                    }
                    try {
                        cell.decode_conserved(0, 1, 0.0);
                    }
                    catch (FlowSolverException e) {
                        writefln("Failed attempt %d: dt= %e", attempt, dt);
                        failedAttempt = true;
                        dt = 0.1*dt;
                        break;
                    }
                    cellCount += nConserved;
                }
            }

            if ( failedAttempt )
                continue;

            // If we get here, things are good. Put flow state into U[0]
            // ready for next iteration.
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (cell; blk.cells) {
                    swap(cell.U[0], cell.U[1]);
                }
            }
            // If we got here, we can break the attempts loop.
            break;
        }
        if ( failedAttempt ) {
            writefln("Step failed: %d", step);
            writeln("Bailing out!");
            exit(1);
        }

        pseudoSimTime += dt;
        wallClockElapsed = 1.0e-3*(Clock.currTime() - wallClockStart).total!"msecs"();  

        // Check on some stopping criteria
        if ( step == nsteps ) {
            writeln("STOPPING: Reached maximum number of steps.");
            finalStep = true;
        }
        if ( normNew <= absGlobalResidReduction && step > nStartUpSteps ) {
            writeln("STOPPING: The absolute global residual is below target value.");
            writefln("          current value= %.12e   target value= %.12e", normNew, absGlobalResidReduction);
            finalStep = true;
        }
        if ( normNew/normRef <= relGlobalResidReduction && step > nStartUpSteps ) {
            writeln("STOPPING: The relative global residual is below target value.");
            writefln("          current value= %.12e   target value= %.12e", normNew/normRef, relGlobalResidReduction);
            finalStep = true;
        }
        if (GlobalConfig.halt_now == 1) {
            writeln("STOPPING: Halt set in control file.");
            finalStep = true;
        }
        
        // Now do some output and diagnostics work
        if ( (step % writeDiagnosticsCount) == 0 || finalStep ) {
            cfl = determine_min_cfl(dt);
            // Write out residuals
            if ( !residualsUpToDate ) {
                max_residuals(currResiduals);
                residualsUpToDate = true;
            }
            fResid = File(residFname, "a");
            fResid.writef("%8d  %20.16e  %20.16e %20.16e %20.16e %3d %5d %.8f %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  %20.16e  ",
                          step, pseudoSimTime, dt, cfl, eta, nRestarts, fnCount, wallClockElapsed, 
                          normNew, normNew/normRef,
                          currResiduals.mass, currResiduals.mass/maxResiduals.mass,
                          currResiduals.momentum.x, currResiduals.momentum.x/maxResiduals.momentum.x,
                          currResiduals.momentum.y, currResiduals.momentum.y/maxResiduals.momentum.y);
            if ( GlobalConfig.dimensions == 3 )
                fResid.writef("%20.16e  %20.16e  ", currResiduals.momentum.z, currResiduals.momentum.z/maxResiduals.momentum.z);
            fResid.writef("%20.16e  %20.16e  ",
                          currResiduals.total_energy, currResiduals.total_energy/maxResiduals.total_energy);
            if ( with_k_omega ) {
                fResid.writef("%20.16e  %20.16e  %20.16e  %20.16e  ",
                              currResiduals.tke, currResiduals.tke/maxResiduals.tke,
                              currResiduals.omega, currResiduals.omega/maxResiduals.omega);
            }
            fResid.write("\n");
            fResid.close();
        }

        // write out the loads
        if ( (step % writeLoadsCount) == 0 || finalStep ) {
            write_boundary_loads_to_file(pseudoSimTime, step);
        }
        
        if ( (step % GlobalConfig.print_count) == 0 || finalStep ) {
            cfl = determine_min_cfl(dt);
            if ( !residualsUpToDate ) {
                max_residuals(currResiduals);
                residualsUpToDate = true;
            }
            auto writer = appender!string();

            formattedWrite(writer, "STEP= %7d  pseudo-time=%10.3e dt=%10.3e cfl=%10.3e  WC=%.1f \n", step, pseudoSimTime, dt, cfl, wallClockElapsed);
            formattedWrite(writer, "RESIDUALS        absolute        relative\n");
            formattedWrite(writer, "  global         %10.6e    %10.6e\n", normNew, normNew/normRef);
            formattedWrite(writer, "  mass           %10.6e    %10.6e\n", currResiduals.mass, currResiduals.mass/maxResiduals.mass);
            formattedWrite(writer, "  x-mom          %10.6e    %10.6e\n", currResiduals.momentum.x, currResiduals.momentum.x/maxResiduals.momentum.x);
            formattedWrite(writer, "  y-mom          %10.6e    %10.6e\n", currResiduals.momentum.y, currResiduals.momentum.y/maxResiduals.momentum.y);
            if ( GlobalConfig.dimensions == 3 )
                formattedWrite(writer, "  z-mom          %10.6e    %10.6e\n", currResiduals.momentum.z, currResiduals.momentum.z/maxResiduals.momentum.z);
            formattedWrite(writer, "  total-energy   %10.6e    %10.6e\n", currResiduals.total_energy, currResiduals.total_energy/maxResiduals.total_energy);
            if ( with_k_omega ) {
                formattedWrite(writer, "  tke            %10.6e    %10.6e\n", currResiduals.tke, currResiduals.tke/maxResiduals.tke);
                formattedWrite(writer, "  omega          %10.6e    %10.6e\n", currResiduals.omega, currResiduals.omega/maxResiduals.omega);
            }
            writeln(writer.data);
        }

        // Write out the flow field, if required
        if ( (step % snapshotsCount) == 0 || finalStep ) {
            if ( !residualsUpToDate ) {
                max_residuals(currResiduals);
                residualsUpToDate = true;
            }
            writefln("-----------------------------------------------------------------------");
            writefln("Writing flow solution at step= %4d; pseudo-time= %6.3e", step, pseudoSimTime);
            writefln("-----------------------------------------------------------------------\n");
            nWrittenSnapshots++;
            if ( nWrittenSnapshots <= nTotalSnapshots ) {
                ensure_directory_is_present(make_path_name!"flow"(nWrittenSnapshots));
                foreach ( iblk, blk; localFluidBlocks ) {
                    // [TODO] PJ 2017-09-02 need to set the file extension properly (here and a few lines below).
                    auto fileName = make_file_name!"flow"(jobName, to!int(iblk), nWrittenSnapshots, "gz");
                    blk.write_solution(fileName, pseudoSimTime);
                }
                restartInfo.pseudoSimTime = pseudoSimTime;
                restartInfo.dt = dt;
                restartInfo.step = step;
                restartInfo.globalResidual = normNew;
                restartInfo.residuals = currResiduals;
                times ~= restartInfo;
                rewrite_times_file(times);
            }
            else {
                // We need to shuffle all of the snapshots...
                foreach ( iSnap; 2 .. nTotalSnapshots+1) {
                    foreach ( iblk; 0 .. localFluidBlocks.length ) {
                        auto fromName = make_file_name!"flow"(jobName, to!int(iblk), iSnap, "gz");
                        auto toName = make_file_name!"flow"(jobName, to!int(iblk), iSnap-1, "gz");
                        rename(fromName, toName);
                    }
                }
                // ... and add the new snapshot.
                foreach ( iblk, blk; localFluidBlocks ) {
                    auto fileName = make_file_name!"flow"(jobName, to!int(iblk), nTotalSnapshots, "gz");        
                    blk.write_solution(fileName, pseudoSimTime);
                }
                remove(times, 1);
                restartInfo.pseudoSimTime = pseudoSimTime;
                restartInfo.dt = dt;
                restartInfo.step = step;
                restartInfo.globalResidual = normNew;
                restartInfo.residuals = currResiduals;
                times[$-1] = restartInfo;
                rewrite_times_file(times);
            }
        }

        if ( finalStep ) break;
        
        if ( !inexactNewtonPhase && normNew/normRef < tau ) {
            // Switch to inexactNewtonPhase
            inexactNewtonPhase = true;
        }

        // Choose a new timestep and eta value.
        auto normRatio = normOld/normNew;
        if ( inexactNewtonPhase ) {
            if ( step < nStartUpSteps ) {
                // Let's assume we're still letting the shock settle
                // when doing low order steps, so we use a power of 0.75
                dtTrial = dt*pow(normOld/normNew, 0.75);
            }
            else {
                // We use a power of 1.0
                dtTrial = dt*normRatio;
            }
            // Apply safeguards to dt
            dtTrial = min(dtTrial, 2.0*dt);
            dtTrial = max(dtTrial, 0.1*dt);
            dt = dtTrial;
        }

        if ( step == nStartUpSteps ) {
            // At the swap-over point from start-up phase to main phase
            // we need to do a few special things.
            // 1. Reset dt to user's choice for this new phase based on cfl1.
            writefln("step= %d dt= %e  cfl1= %f", step, dt, cfl1);
            dt = determine_initial_dt(cfl1);
            writefln("after choosing new timestep: %e", dt);
            // 2. Reset the inexact Newton phase.
            //    We'll take some constant timesteps at the new dt
            //    until the residuals have dropped.
            inexactNewtonPhase = false;
        }

        if ( step > nStartUpSteps ) {
            // Adjust eta1 according to eta update strategy
            final switch (etaStrategy) {
            case EtaStrategy.constant:
                break;
            case EtaStrategy.geometric:
                eta1 *= etaRatioPerStep;
                // but don't go past eta1_min
                eta1 = fmax(eta1, eta1_min);
                break;
            case EtaStrategy.adaptive, EtaStrategy.adaptive_capped:
                // Use Eisenstat & Walker adaptive strategy no. 2
                auto etaOld = eta1;
                eta1 = gamma*pow(1./normRatio, alpha);
                // Apply Eisenstat & Walker safeguards
                auto testVal = gamma*pow(etaOld, alpha);
                if ( testVal > 0.1 ) {
                    eta1 = fmax(eta1, testVal);
                }
                // and don't let eta1 get larger than a max value
                eta1 = fmin(eta1, eta1_max);
                if ( etaStrategy == EtaStrategy.adaptive_capped ) {
                    // Additionally, we cap the maximum so that we never
                    // retreat on the eta1 value.
                    eta1_max = fmin(eta1, eta1_max);
                    // but don't let our eta1_max cap get too tiny
                    eta1_max = fmax(eta1, eta1_min);
                }
                break;
            }
        }

        normOld = normNew;
    }

}

void allocate_global_workspace()
{
    size_t mOuter = to!size_t(GlobalConfig.sssOptions.maxOuterIterations);
    g0.length = mOuter+1;
    g1.length = mOuter+1;
    h.length = mOuter+1;
    hR.length = mOuter+1;
    H0 = new Matrix!number(mOuter+1, mOuter);
    H1 = new Matrix!number(mOuter+1, mOuter);
    Gamma = new Matrix!number(mOuter+1, mOuter+1);
    Q0 = new Matrix!number(mOuter+1, mOuter+1);
    Q1 = new Matrix!number(mOuter+1, mOuter+1);
}

double determine_initial_dt(double cflInit)
{
    double signal, dt_local, dt;
    bool first = true;

    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            signal = cell.signal_frequency();
            dt_local = cflInit / signal;
            if ( first )
                dt = dt_local;
            else
                dt = fmin(dt, dt_local);
        }
    }
    return dt;
}

double determine_min_cfl(double dt)
{
    double signal, cfl_local, cfl;
    bool first = true;

    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) {
            signal = cell.signal_frequency();
            cfl_local = dt * signal;
            if ( first )
                cfl = cfl_local;
            else
                cfl = fmin(cfl, cfl_local);
        }
    }
    return cfl;
}

void evalRHS(double pseudoSimTime, int ftl)
{
    fnCount++;
    
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
    }
    exchange_ghost_cell_boundary_data(pseudoSimTime, 0, ftl);
    foreach (blk; localFluidBlocks) {
        blk.applyPreReconAction(pseudoSimTime, 0, ftl);
    }
    
    // We don't want to switch between flux calculator application while
    // doing the Frechet derivative, so we'll only search for shock points
    // at ftl = 0, which is when the F(U) evaluation is made.
    if ( ftl == 0 &&
         (GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv || FluxCalculator.adaptive_hlle_ausmdv) ) {
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.detect_shock_points();
        }
    }
    
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase0();
    }
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.convective_flux_phase1();
    }
    foreach (blk; localFluidBlocks) {
        blk.applyPostConvFluxAction(pseudoSimTime, 0, ftl);
    }
    if (GlobalConfig.viscous) {
        foreach (blk; localFluidBlocks) {
            blk.applyPreSpatialDerivActionAtBndryFaces(pseudoSimTime, 0, ftl);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.flow_property_spatial_derivatives(0); 
            blk.estimate_turbulence_viscosity();
            blk.viscous_flux();
        }
        foreach (blk; localFluidBlocks) {
            blk.applyPostDiffFluxAction(pseudoSimTime, 0, ftl);
        }
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        foreach (i, cell; blk.cells) {
            cell.add_inviscid_source_vector(0, 0.0);
            if (blk.myConfig.viscous) {
                cell.add_viscous_source_vector(local_with_k_omega);
            }
            if (blk.myConfig.udf_source_terms) {
                addUDFSourceTermsToCell(blk.myL, cell, 0, 
                                        pseudoSimTime, blk.myConfig.gmodel);
            }
            cell.time_derivatives(0, ftl, local_with_k_omega);
        }
    }
}

void evalJacobianVecProd(double pseudoSimTime, double sigma)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;
    size_t OMEGA = omegaIdx;

    // We perform a Frechet derivative to evaluate J*D^(-1)v
    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        blk.clear_fluxes_of_conserved_quantities();
        foreach (cell; blk.cells) cell.clear_source_vector();
        int cellCount = 0;
        foreach (cell; blk.cells) {
            cell.U[1].copy_values_from(cell.U[0]);
            cell.U[1].mass += sigma*blk.maxRate.mass*blk.zed[cellCount+MASS];
            cell.U[1].momentum.refx += sigma*blk.maxRate.momentum.x*blk.zed[cellCount+X_MOM];
            cell.U[1].momentum.refy += sigma*blk.maxRate.momentum.y*blk.zed[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                cell.U[1].momentum.refz += sigma*blk.maxRate.momentum.z*blk.zed[cellCount+Z_MOM];
            cell.U[1].total_energy += sigma*blk.maxRate.total_energy*blk.zed[cellCount+TOT_ENERGY];
            if ( local_with_k_omega ) {
                cell.U[1].tke += sigma*blk.maxRate.tke*blk.zed[cellCount+TKE];
                cell.U[1].omega += sigma*blk.maxRate.omega*blk.zed[cellCount+OMEGA];
            }
            cell.decode_conserved(0, 1, 0.0);
            cellCount += nConserved;
        }
    }
    evalRHS(pseudoSimTime, 1);
    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.zed[cellCount+MASS] = (-cell.dUdt[1].mass - blk.FU[cellCount+MASS])/(sigma*blk.maxRate.mass);
            blk.zed[cellCount+X_MOM] = (-cell.dUdt[1].momentum.x - blk.FU[cellCount+X_MOM])/(sigma*blk.maxRate.momentum.x);
            blk.zed[cellCount+Y_MOM] = (-cell.dUdt[1].momentum.y - blk.FU[cellCount+Y_MOM])/(sigma*blk.maxRate.momentum.y);
            if ( blk.myConfig.dimensions == 3 )
                blk.zed[cellCount+Z_MOM] = (-cell.dUdt[1].momentum.z - blk.FU[cellCount+Z_MOM])/(sigma*blk.maxRate.momentum.z);
            blk.zed[cellCount+TOT_ENERGY] = (-cell.dUdt[1].total_energy - blk.FU[cellCount+TOT_ENERGY])/(sigma*blk.maxRate.total_energy);
            if ( local_with_k_omega ) {
                blk.zed[cellCount+TKE] = (-cell.dUdt[1].tke - blk.FU[cellCount+TKE])/(sigma*blk.maxRate.tke);
                blk.zed[cellCount+OMEGA] = (-cell.dUdt[1].omega - blk.FU[cellCount+OMEGA])/(sigma*blk.maxRate.omega);
            }
            cellCount += nConserved;
        }
    }
}


string dot_over_blocks(string dot, string A, string B)
{
    return `
foreach (blk; parallel(localFluidBlocks,1)) {
   blk.dotAcc = 0.0;
   foreach (k; 0 .. blk.nvars) {
      blk.dotAcc += blk.`~A~`[k].re*blk.`~B~`[k].re;
   }
}
`~dot~` = 0.0;
foreach (blk; localFluidBlocks) `~dot~` += blk.dotAcc;`;

}

string norm2_over_blocks(string norm2, string blkMember)
{
    return `
foreach (blk; parallel(localFluidBlocks,1)) {
   blk.normAcc = 0.0;
   foreach (k; 0 .. blk.nvars) {
      blk.normAcc += blk.`~blkMember~`[k].re*blk.`~blkMember~`[k].re;
   }
}
`~norm2~` = 0.0;
foreach (blk; localFluidBlocks) `~norm2~` += blk.normAcc;
`~norm2~` = sqrt(`~norm2~`);`;

}

void rpcGMRES_solve(double pseudoSimTime, double dt, double eta, double sigma, bool usePreconditioner,
                    ref double residual, ref int nRestarts)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;
    size_t MASS = massIdx;
    size_t X_MOM = xMomIdx;
    size_t Y_MOM = yMomIdx;
    size_t Z_MOM = zMomIdx;
    size_t TOT_ENERGY = totEnergyIdx;
    size_t TKE = tkeIdx;
    size_t OMEGA = omegaIdx;

    number resid;
    
    int interpOrderSave = GlobalConfig.interpolation_order;
    // Presently, just do one block
    int maxIters = GlobalConfig.sssOptions.maxOuterIterations;
    // We add 1 because the user thinks of "re"starts, so they
    // might legitimately ask for no restarts. We still have
    // to execute at least once.
    int maxRestarts = GlobalConfig.sssOptions.maxRestarts + 1; 
    size_t m = to!size_t(maxIters);
    size_t r;
    size_t iterCount;

    // Variables for max rates of change
    // Use these for equation scaling.
    double minNonDimVal = 1.0; // minimum value used for non-dimensionalisation
                               // when our time rates of change are very small
                               // then we'll avoid non-dimensionalising by 
                               // values close to zero.

    // 1. Evaluate r0, beta, v1
    evalRHS(pseudoSimTime, 0);
    // Store dUdt[0] as F(U)
    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        int cellCount = 0;
        blk.maxRate.mass = 0.0;
        blk.maxRate.momentum.refx = 0.0;
        blk.maxRate.momentum.refy = 0.0;
        if ( blk.myConfig.dimensions == 3 )
            blk.maxRate.momentum.refz = 0.0;
        blk.maxRate.total_energy = 0.0;
        if ( local_with_k_omega ) {
            blk.maxRate.tke = 0.0;
            blk.maxRate.omega = 0.0;
        }
        foreach (cell; blk.cells) {
            blk.FU[cellCount+MASS] = -cell.dUdt[0].mass;
            blk.FU[cellCount+X_MOM] = -cell.dUdt[0].momentum.x;
            blk.FU[cellCount+Y_MOM] = -cell.dUdt[0].momentum.y;
            if ( blk.myConfig.dimensions == 3 )
                blk.FU[cellCount+Z_MOM] = -cell.dUdt[0].momentum.z;
            blk.FU[cellCount+TOT_ENERGY] = -cell.dUdt[0].total_energy;
            if ( local_with_k_omega ) {
                blk.FU[cellCount+TKE] = -cell.dUdt[0].tke;
                blk.FU[cellCount+OMEGA] = -cell.dUdt[0].omega;
            }
            cellCount += nConserved;

            blk.maxRate.mass = fmax(blk.maxRate.mass, fabs(cell.dUdt[0].mass));
            blk.maxRate.momentum.refx = fmax(blk.maxRate.momentum.x, fabs(cell.dUdt[0].momentum.x));
            blk.maxRate.momentum.refy = fmax(blk.maxRate.momentum.y, fabs(cell.dUdt[0].momentum.y));
            if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.momentum.refz = fmax(blk.maxRate.momentum.z, fabs(cell.dUdt[0].momentum.z));
            blk.maxRate.total_energy = fmax(blk.maxRate.total_energy, fabs(cell.dUdt[0].total_energy));
            if ( local_with_k_omega ) {
                blk.maxRate.tke = fmax(blk.maxRate.tke, fabs(cell.dUdt[0].tke));
                blk.maxRate.omega = fmax(blk.maxRate.omega, fabs(cell.dUdt[0].omega));
            }
        }
    }
    number maxMass = 0.0;
    number maxMomX = 0.0;
    number maxMomY = 0.0;
    number maxMomZ = 0.0;
    number maxEnergy = 0.0;
    number maxTke = 0.0;
    number maxOmega = 0.0;
    foreach (blk; localFluidBlocks) {
        maxMass = fmax(maxMass, blk.maxRate.mass);
        maxMomX = fmax(maxMomX, blk.maxRate.momentum.x);
        maxMomY = fmax(maxMomY, blk.maxRate.momentum.y);
        if ( blk.myConfig.dimensions == 3 )
            maxMomZ = fmax(maxMomZ, blk.maxRate.momentum.z);
        maxEnergy = fmax(maxEnergy, blk.maxRate.total_energy);
        if ( with_k_omega ) {
            maxTke = fmax(maxTke, blk.maxRate.tke);
            maxOmega = fmax(maxOmega, blk.maxRate.omega);
        }
    }
    // Place some guards when time-rate-of-changes are very small.
    maxMass = fmax(maxMass, minNonDimVal);
    maxMomX = fmax(maxMomX, minNonDimVal);
    maxMomY = fmax(maxMomY, minNonDimVal);
    if ( GlobalConfig.dimensions == 3 )
        maxMomZ = fmax(maxMomZ, minNonDimVal);
    maxEnergy = fmax(maxEnergy, minNonDimVal);
    number maxMom = max(maxMomX, maxMomY, maxMomZ);
    if ( with_k_omega ) {
        maxTke = fmax(maxTke, minNonDimVal);
        maxOmega = fmax(maxOmega, minNonDimVal);
    }

    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.maxRate.mass = maxMass;
        blk.maxRate.momentum.refx = maxMomX;
        blk.maxRate.momentum.refy = maxMomY;
        if ( blk.myConfig.dimensions == 3 )
                blk.maxRate.momentum.refz = maxMomZ;
        blk.maxRate.total_energy = maxEnergy;
        if ( with_k_omega ) {
            blk.maxRate.tke = maxTke;
            blk.maxRate.omega = maxOmega;
        }
    }
    double unscaledNorm2;
    mixin(norm2_over_blocks("unscaledNorm2", "FU"));

    // Initialise some arrays and matrices that have already been allocated
    g0[] = to!number(0.0);
    g1[] = to!number(0.0);
    H0.zeros();
    H1.zeros();

    double dtInv = 1.0/dt;

    // We'll scale r0 against these max rates of change.
    // r0 = b - A*x0
    // Taking x0 = [0] (as is common) gives r0 = b = FU
    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        blk.x0[] = to!number(0.0);
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.r0[cellCount+MASS] = -(1./blk.maxRate.mass)*blk.FU[cellCount+MASS];
            blk.r0[cellCount+X_MOM] = -(1./blk.maxRate.momentum.x)*blk.FU[cellCount+X_MOM];
            blk.r0[cellCount+Y_MOM] = -(1./blk.maxRate.momentum.y)*blk.FU[cellCount+Y_MOM];
            if ( blk.myConfig.dimensions == 3 )
                blk.r0[cellCount+Z_MOM] = -(1./blk.maxRate.momentum.z)*blk.FU[cellCount+Z_MOM];
            blk.r0[cellCount+TOT_ENERGY] = -(1./blk.maxRate.total_energy)*blk.FU[cellCount+TOT_ENERGY];
            if ( local_with_k_omega ) {
                blk.r0[cellCount+TKE] = -(1./blk.maxRate.tke)*blk.FU[cellCount+TKE];
                blk.r0[cellCount+OMEGA] = -(1./blk.maxRate.omega)*blk.FU[cellCount+OMEGA];
            }
            cellCount += nConserved;
        }
    }
    // Then compute v = r0/||r0||
    number betaTmp;
    mixin(norm2_over_blocks("betaTmp", "r0"));
    number beta = betaTmp;
    g0[0] = beta;
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (k; 0 .. blk.nvars) {
            blk.v[k] = blk.r0[k]/beta;
            blk.V[k,0] = blk.v[k];
        }
    }

    // Compute tolerance
    auto outerTol = eta*beta;

    // 2. Start outer-loop of restarted GMRES
    for ( r = 0; r < maxRestarts; r++ ) {
        // 2a. Begin iterations
        foreach (j; 0 .. m) {
            iterCount = j+1;
            if (usePreconditioner) {
                bool savedViscousFlag = GlobalConfig.viscous;
                GlobalConfig.viscous = false;
                foreach (blk; parallel(localFluidBlocks,1)) {
                    if (r == 0 && j == 0) {
                        blk.myConfig.viscous = false;
                        sss_preconditioner(blk, nConserved, dt, 1.0e-6, 1.0e-6);
                        blk.myConfig.viscous = savedViscousFlag;
                    }
                    int cellCount = 0;
                    number[] tmp;
                    tmp.length = nConserved;
                    foreach (cell; blk.cells) {
                        LUSolve!number(cell.dConservative, cell.pivot,
                                       blk.v[cellCount..cellCount+nConserved], tmp);
                        blk.zed[cellCount..cellCount+nConserved] = tmp[];
                        cellCount += nConserved;
                    }
                }
                GlobalConfig.viscous = savedViscousFlag;
            }
            else {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    blk.zed[] = blk.v[];
                }
            }
            // Prepare 'w' with (I/dt)(P^-1)v term;
            foreach (blk; parallel(localFluidBlocks,1)) {
                double dtInv = 1.0/dt;
                foreach (idx; 0..blk.w.length) blk.w[idx] = dtInv*blk.zed[idx];
            }
            
            // Evaluate Jz and place in z
            evalJacobianVecProd(pseudoSimTime, sigma);

            // Now we can complete calculation of w
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars)  blk.w[k] += blk.zed[k];
            }
            // The remainder of the algorithm looks a lot like any standard
            // GMRES implementation (for example, see smla.d)
            foreach (i; 0 .. j+1) {
                foreach (blk; parallel(localFluidBlocks,1)) {
                    // Extract column 'i'
                    foreach (k; 0 .. blk.nvars ) blk.v[k] = blk.V[k,i]; 
                }
                number H0_ij_tmp;
                mixin(dot_over_blocks("H0_ij_tmp", "w", "v"));
                number H0_ij = H0_ij_tmp;
                H0[i,j] = H0_ij;
                foreach (blk; parallel(localFluidBlocks,1)) {
                    foreach (k; 0 .. blk.nvars) blk.w[k] -= H0_ij*blk.v[k]; 
                }
            }
            number H0_jp1j_tmp;
            mixin(norm2_over_blocks("H0_jp1j_tmp", "w"));
            number H0_jp1j = H0_jp1j_tmp;
            H0[j+1,j] = H0_jp1j;
        
            foreach (blk; parallel(localFluidBlocks,1)) {
                foreach (k; 0 .. blk.nvars) {
                    blk.v[k] = blk.w[k]/H0_jp1j;
                    blk.V[k,j+1] = blk.v[k];
                }
            }

            // Build rotated Hessenberg progressively
            if ( j != 0 ) {
                // Extract final column in H
                foreach (i; 0 .. j+1) h[i] = H0[i,j];
                // Rotate column by previous rotations (stored in Q0)
                nm.bbla.dot(Q0, j+1, j+1, h, hR);
                // Place column back in H
                foreach (i; 0 .. j+1) H0[i,j] = hR[i];
            }
            // Now form new Gamma
            Gamma.eye();
            auto denom = sqrt(H0[j,j]*H0[j,j] + H0[j+1,j]*H0[j+1,j]);
            auto s_j = H0[j+1,j]/denom; 
            auto c_j = H0[j,j]/denom;
            Gamma[j,j] = c_j; Gamma[j,j+1] = s_j;
            Gamma[j+1,j] = -s_j; Gamma[j+1,j+1] = c_j;
            // Apply rotations
            nm.bbla.dot(Gamma, j+2, j+2, H0, j+1, H1);
            nm.bbla.dot(Gamma, j+2, j+2, g0, g1);
            // Accumulate Gamma rotations in Q.
            if ( j == 0 ) {
                copy(Gamma, Q1);
            }
            else {
                nm.bbla.dot!number(Gamma, j+2, j+2, Q0, j+2, Q1);
            }

            // Prepare for next step
            copy(H1, H0);
            g0[] = g1[];
            copy(Q1, Q0);

            // Get residual
            resid = fabs(g1[j+1]);
            // DEBUG:
            //      writefln("OUTER: restart-count= %d iteration= %d, resid= %e", r, j, resid);
            if ( resid <= outerTol ) {
                m = j+1;
                // DEBUG:
                //      writefln("OUTER: TOL ACHIEVED restart-count= %d iteration-count= %d, resid= %e", r, m, resid);
                break;
            }
        }

        if (iterCount == maxIters)
            m = maxIters;

        // At end H := R up to row m
        //        g := gm up to row m
        upperSolve!number(H1, to!int(m), g1);
        // In serial, distribute a copy of g1 to each block
        foreach (blk; localFluidBlocks) blk.g1[] = g1[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot!number(blk.V, blk.nvars, m, blk.g1, blk.zed);
        }
        
        if (usePreconditioner) {
            foreach(blk; parallel(localFluidBlocks,1)) {
                int cellCount = 0;
                number[] tmp;
                tmp.length = nConserved;
                foreach (cell; blk.cells) {
                    LUSolve!number(cell.dConservative, cell.pivot,
                                   blk.zed[cellCount..cellCount+nConserved], tmp);
                    blk.dU[cellCount..cellCount+nConserved] = tmp[];
                    cellCount += nConserved;
                }
            }
        }
        else {
            foreach(blk; parallel(localFluidBlocks,1)) {
                blk.dU[] = blk.zed[];
            }
        }
        
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) blk.dU[k] += blk.x0[k];
        }

        if ( resid <= outerTol || r+1 == maxRestarts ) {
            // DEBUG:  writefln("resid= %e outerTol= %e  r+1= %d  maxRestarts= %d", resid, outerTol, r+1, maxRestarts);
            // DEBUG:  writefln("Breaking restart loop.");
            break;
        }

        // Else, we prepare for restart by setting x0 and computing r0
        // Computation of r0 as per Fraysee etal (2005)
        foreach (blk; parallel(localFluidBlocks,1)) {
            blk.x0[] = blk.dU[];
        }

        foreach (blk; localFluidBlocks) copy(Q1, blk.Q1);
        // Set all values in g0 to 0.0 except for final (m+1) value
        foreach (i; 0 .. m) g0[i] = 0.0;
        foreach (blk; localFluidBlocks) blk.g0[] = g0[];
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot(blk.Q1, m, m+1, blk.g0, blk.g1);
        }
        foreach (blk; parallel(localFluidBlocks,1)) {
            nm.bbla.dot(blk.V, blk.nvars, m+1, blk.g1, blk.r0);
        }
        
        mixin(norm2_over_blocks("betaTmp", "r0"));
        beta = betaTmp;
        // DEBUG: writefln("OUTER: ON RESTART beta= %e", beta);
        foreach (blk; parallel(localFluidBlocks,1)) {
            foreach (k; 0 .. blk.nvars) {
                blk.v[k] = blk.r0[k]/beta;
                blk.V[k,0] = blk.v[k];
            }
        }
        // Re-initialise some vectors and matrices for restart
        g0[] = to!number(0.0);
        g1[] = to!number(0.0);
        H0.zeros();
        H1.zeros();
        // And set first residual entry
        g0[0] = beta;

    }
    // Rescale
    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        int cellCount = 0;
        foreach (cell; blk.cells) {
            blk.dU[cellCount+MASS] *= blk.maxRate.mass;
            blk.dU[cellCount+X_MOM] *= blk.maxRate.momentum.x;
            blk.dU[cellCount+Y_MOM] *= blk.maxRate.momentum.y;
            if ( blk.myConfig.dimensions == 3 )
                blk.dU[cellCount+Z_MOM] *= blk.maxRate.momentum.z;
            blk.dU[cellCount+TOT_ENERGY] *= blk.maxRate.total_energy;
            if ( local_with_k_omega ) {
                blk.dU[cellCount+TKE] *= blk.maxRate.tke;
                blk.dU[cellCount+OMEGA] *= blk.maxRate.omega;
            }
            cellCount += nConserved;
        }
    }

    residual = unscaledNorm2;
    nRestarts = to!int(r);
}

void max_residuals(ConservedQuantities residuals)
{
    // Make a stack-local copy of conserved quantities info
    size_t nConserved = nConservedQuantities;

    foreach (blk; parallel(localFluidBlocks,1)) {
        bool local_with_k_omega = with_k_omega;
        blk.residuals.copy_values_from(blk.cells[0].dUdt[0]);
        blk.residuals.mass = fabs(blk.residuals.mass);
        blk.residuals.momentum.refx = fabs(blk.residuals.momentum.x);
        blk.residuals.momentum.refy = fabs(blk.residuals.momentum.y);
        if ( blk.myConfig.dimensions == 3 )
            blk.residuals.momentum.refz = fabs(blk.residuals.momentum.z);
        blk.residuals.total_energy = fabs(blk.residuals.total_energy);
        if ( local_with_k_omega ) {
            blk.residuals.tke = fabs(blk.residuals.tke);
            blk.residuals.omega = fabs(blk.residuals.omega);
        }
        number massLocal, xMomLocal, yMomLocal, zMomLocal, energyLocal, tkeLocal, omegaLocal;
        foreach (cell; blk.cells) {
            massLocal = cell.dUdt[0].mass;
            xMomLocal = cell.dUdt[0].momentum.x;
            yMomLocal = cell.dUdt[0].momentum.y;
            zMomLocal = cell.dUdt[0].momentum.z;
            energyLocal = cell.dUdt[0].total_energy;
            if ( local_with_k_omega ) {
                tkeLocal = cell.dUdt[0].tke;
                omegaLocal = cell.dUdt[0].omega;
            }
            blk.residuals.mass = fmax(blk.residuals.mass, massLocal);
            blk.residuals.momentum.refx = fmax(blk.residuals.momentum.x, xMomLocal);
            blk.residuals.momentum.refy = fmax(blk.residuals.momentum.y, yMomLocal);
            if ( blk.myConfig.dimensions == 3 )
                blk.residuals.momentum.refz = fmax(blk.residuals.momentum.z, zMomLocal);
            blk.residuals.total_energy = fmax(blk.residuals.total_energy, energyLocal);
            if ( local_with_k_omega ) {
                blk.residuals.tke = fmax(blk.residuals.tke, tkeLocal);
                blk.residuals.omega = fmax(blk.residuals.omega, omegaLocal);
            }
        }
    }
    residuals.copy_values_from(localFluidBlocks[0].residuals);
    foreach (blk; localFluidBlocks) {
        residuals.mass = fmax(residuals.mass, blk.residuals.mass);
        residuals.momentum.refx = fmax(residuals.momentum.x, blk.residuals.momentum.x);
        residuals.momentum.refy = fmax(residuals.momentum.y, blk.residuals.momentum.y);
        if ( blk.myConfig.dimensions == 3 )
            residuals.momentum.refz = fmax(residuals.momentum.z, blk.residuals.momentum.z);
        residuals.total_energy = fmax(residuals.total_energy, blk.residuals.total_energy);
        if ( with_k_omega ) {
            residuals.tke = fmax(residuals.tke, blk.residuals.tke);
            residuals.omega = fmax(residuals.omega, blk.residuals.omega);
        }
    }
}

void rewrite_times_file(RestartInfo[] times)
{
    auto fname = format("%s.times", GlobalConfig.base_file_name);
    auto f = File(fname, "w");
    f.writeln("# tindx sim_time dt_global");
    foreach (i, rInfo; times) {
        if ( GlobalConfig.dimensions == 2 ) {
            f.writef("%04d %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e",
                     i, rInfo.pseudoSimTime, rInfo.dt, rInfo.step,
                     rInfo.globalResidual, rInfo.residuals.mass,
                     rInfo.residuals.momentum.x, rInfo.residuals.momentum.y,
                     rInfo.residuals.total_energy);
        }
        else {
            f.writef("%04d %.18e %.18e %d %.18e %.18e %.18e %.18e %.18e %.18e",
                     i, rInfo.pseudoSimTime, rInfo.dt, rInfo.step,
                     rInfo.globalResidual, rInfo.residuals.mass,
                     rInfo.residuals.momentum.x, rInfo.residuals.momentum.y, rInfo.residuals.momentum.z,
                     rInfo.residuals.total_energy);
        }
        if ( GlobalConfig.turbulence_model == TurbulenceModel.k_omega ) {
            f.writef(" %.18e %.18e", rInfo.residuals.tke, rInfo.residuals.omega);
        }
        f.write("\n");
    }
    f.close();
}
