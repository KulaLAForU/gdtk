// characteristic_bc.d
// 
// JSW 2026-08-25
// 
// Idea came from Christian Stemmer
// Formulation of characteristic boundary condition at the free stream,
// Original from Paul Harris' thesis UofA
// Only working for structured mesh

module bc.ghost_cell_effect.characteristic_bc;

import std.conv;
import std.json;
import std.math;
import std.stdio;
import std.string;

import gas;
import geom;

import lmr.bc;
import lmr.flowstate;
import lmr.fluidblock;
import lmr.fluidfvcell;
import lmr.fvinterface;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.sfluidblock;

class GhostCellCharacteristic : GhostCellEffect {
public:
    // Perhaps here need distinguish this BC is for north or east by:
    // BoundaryCondition bc = blk.bc[which_boundary];
    //     if (bc.outsigns[f.i_bndry] == 1) and ()
    @nogc
    override void apply_apply_structured_grid(double t, int gtl, int ftl)
    {
        auto blk = cast(SFluidBlock) this.blk;
        assert(blk !is null, "Oops, this should be an SFluidBlock object.");
        assert(blk.n_ghost_cell_layers == 2, "Oops, the ghost cell layers should be 2 for this characteristic BC");
        BoundaryCondition bc = blk.bc[which_boundary];
        foreach (i, f; bc.faces) {
            foreach (n; 0 .. blk.n_ghost_cell_layers) {
                FluidFVCell Cell1, ghost1, ghost2;
                if (bc.outsigns[i] == 1) {
                    Cell1 = f.left_cells[0]; ghost1 = f.right_cell[0]; ghost2 = f.right_cell[1];
                    characteristic_condition(Cell1.fs, ghost1.fs, ghost2.fs);
                }
            }
        }

    }

private:
    // @nogc
    // number local_mach_number(const FlowState* fs0)
    // {
    //     local_v = sqrt(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y + fs.vel.z*fs.vel.z);
    //     local_m = local_v / fs.gas.a;
    // }
    
    @nogc
    void characteristic_condition(FlowState* fs0, FlowState* fs1, FlowState* fs2)
    {
        foreach (i; 0 .. 2) {
            local_v = sqrt(fs.vel.x*fs.vel.x + fs.vel.y*fs.vel.y + fs.vel.z*fs.vel.z);
            local_m = local_v / fs.gas.a; // local mach number
            fmu =  asin(1/local_m) + atan(fs.vel.v / fs.vel.u); // local mach angle
            // here need pass dx and dy to find the helper points
            dyx = dy/dx;
            fx = dyx / tan(fmu);
            fx1 = n - fx;
            fx2 = n - 2*fx;
            f1 = fx1 - abs(fx1);
            f2 = fx2 - abs(fx2);
            // linear interpolation


        }
    }

    @nogc
    void linear_interpolation()
    {

    }
    
    @nogc
    void characteristic_condition()
    {

    }

}