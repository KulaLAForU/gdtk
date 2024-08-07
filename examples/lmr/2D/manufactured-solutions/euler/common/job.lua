-- mms.lua
--
-- Authors: Rowan G. and Peter J.
-- Date: 2015-03-17
-- History: Ported from eilmer3 example
--          Presently only Euler case on regular grid
--
--          2015-06-03
--          Added scripting to configure case based on 'case.txt'
--
--          2018, KAD
--          Ported to steady-state solver (Kyle A. Damm, 2018)
--
--          2023-07-18, RJG
--          Updated for lmr5, specifically:
--          1. staged prep (using registerFluidGrid); and
--          2. new input for steady-state solver configuration
--

config.dimensions = 2
config.solver_mode = "steady"

-- Case is configured by a higher-level controller.
-- Options are placed in config.txt as simple
-- variable assignment statements and then used
-- in this job script.
dofile("config.txt")

setGasModel('very-viscous-air.lua')
R = 287.0
p0 = 1.0e5; T0 = p0/R;
u0 = 800.0; v0 = 800.0
initial = FlowState:new{p=p0, T=T0, velx=u0, vely=v0}
flowStates = {initial=initial}

p00 = Vector3:new{x=0.0, y=0.0}
p10 = Vector3:new{x=1.0, y=0.0}
p01 = Vector3:new{x=0.0, y=1.0}
p11 = Vector3:new{x=1.0, y=1.0}
nicell = ncells; njcell = ncells
sgrid=StructuredGrid:new{psurface=CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01},
                         niv=nicell+1, njv=njcell+1}
grid = registerFluidGrid{
   grid=UnstructuredGrid:new{sgrid=sgrid},
   fsTag="initial",
   bcTags={[Face.north]="udf", [Face.east]="udf", [Face.south]="udf", [Face.west]="udf"}

}

bcs = {udf=UserDefinedBC:new{fileName='udf-bc.lua'}}
config.apply_bcs_in_parallel = false

makeFluidBlocks(bcs, flowStates)

config.interpolation_order = interpolation_order
config.flux_calculator = flux_calculator
config.udf_source_terms = true
config.udf_source_terms_file = 'udf-source-terms.lua'
config.viscous = false
-- Do NOT use the limiters for the verification tests
config.apply_limiter = false
config.extrema_clipping = false

NewtonKrylovGlobalConfig{
   number_of_steps_for_setting_reference_residuals = 0,
   max_newton_steps = 200,
   stop_on_relative_residual = 1.0e-12,
   number_of_phases = 2,
   max_steps_in_initial_phases = { 100 },
   inviscid_cfl_only = true,
   use_line_search = false,
   use_physicality_check = false,
   max_linear_solver_iterations = 50,
   max_linear_solver_restarts = 0,
   use_scaling = true,
   frechet_derivative_perturbation = 1.0e-50,
   use_preconditioner = true,
   preconditioner_perturbation = 1.0e-50,
   preconditioner = "ilu",
   ilu_fill = 0,
   total_snapshots = 2,
   steps_between_snapshots = 10,
   steps_between_diagnostics = 10
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 1,
   jacobian_interpolation_order = 1,
   frozen_preconditioner = true,
   frozen_limiter_for_jacobian = false,
   use_adaptive_preconditioner = false,
   steps_between_preconditioner_update = 100,
   linear_solve_tolerance = 0.1,
   use_auto_cfl = true,
   threshold_relative_residual_for_cfl_growth = 0.99,
   start_cfl = 1.0,
   max_cfl = 1.0e6,
   auto_cfl_exponent = 1.0
}

NewtonKrylovPhase:new{
   residual_interpolation_order = 2,
   jacobian_interpolation_order = 2,
   start_cfl = 100.0,
}
