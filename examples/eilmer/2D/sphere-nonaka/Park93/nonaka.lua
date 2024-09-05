config.title = "Sphere fired into air."
print(config.title)

no_flow_times = 5.0
Db = 12.7e-3
Rc = Db/2
setGasModel('air-5sp-2T.gas')
p_inf = 666.611842 -- Pa
T_inf = 293.0 -- K
u_inf = 4000.0 -- m/s
inflow = FlowState:new{p=p_inf, T=T_inf, T_modes={T_inf}, velx=u_inf,
                       massf={N2=0.767,O2=0.233}}

body_flow_time = Db/u_inf
config.reacting = true
config.reactions_file = 'air-5sp-6r-2T.chem'
config.energy_exchange_file = 'air-VT.exch'
config.dimensions = 2
config.axisymmetric = true
config.flux_calculator = "ausmdv"
config.interpolation_order = 2
config.interpolation_delay = 3*body_flow_time
-- config.grid_motion = "shock_fitting"
-- config.gasdynamic_update_scheme = "moving_grid_2_stage"
-- config.shock_fitting_delay = body_flow_time
config.gasdynamic_update_scheme = "euler"
config.max_invalid_cells = 20
config.adjust_invalid_cell_data = true
config.max_time = no_flow_times*body_flow_time
config.max_step = 8000000
-- config.dt_init = 1.0e-10
config.dt_init = 1.225e-10
config.cfl_value = 0.5
config.dt_plot = config.max_time/30
config.viscous = true

config.with_local_time_stepping = true

-- config.x_order = 1 -- linear extrapolation in the ghost cells

-- loads
config.write_loads = true
config.boundary_groups_for_loads = "loads"
config.dt_loads = config.dt_plot

a = Vector3:new{x=Rc, y=0.0}
b = Vector3:new{x=0.0, y=0.0}
c = Vector3:new{x=Rc, y=Rc}
d = Vector3:new{x=-0.19*Rc, y=0.0}
e = Vector3:new{x=-0.18*Rc, y=Rc}
f = Vector3:new{x=0.5*Rc,  y=1.4*Rc}
g = Vector3:new{x=Rc,      y=1.75*Rc}
bc = Arc:new{p0=b, p1=c, centre=a}
dg = Bezier:new{points={d, e, f, g}}
db = Line:new{p0=d, p1=b}
gc = Line:new{p0=g, p1=c}

 
psurf = makePatch{north=gc, east=bc, south=db, west=dg}
nx = 100; ny = 150
cf_circum = RobertsFunction:new{end0=false, end1=true, beta=1.3}
-- cf_radial = GeometricFunction:new{a=0.001, r=1.2, N=nx+1, reverse=true}

grid = StructuredGrid:new{psurface=psurf,niv=nx+1, njv=ny+1, cfList={south=cf_circum, north=cf_circum}}
-- grid = StructuredGrid:new{psurface=psurf, niv=nx+1, njv=ny+1, cfList={south=cf_radial, north=cf_radial}}

blk = FBArray:new{grid=grid,initialState=inflow, label='blk',
                  bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
                          east=WallBC_NoSlip_FixedT:new{Twall=1000, group="loads", catalytic_type = "equilibrium" },
                          north=OutFlowBC_Simple:new{}},
                  nib=1, njb=16}



-- local billig_patch = require "billig_patch"
-- M_inf = u_inf/inflow.a
-- print("M_inf=", M_inf)
-- bp = billig_patch.make_patch{Minf=M_inf, R=Rc, scale=1.0}
-- cf_circum = RobertsFunction:new{end0=true, end1=false, beta=1.1}
-- grid = StructuredGrid:new{psurface=bp.patch, niv=61, njv=61,
-- 			  cfList={west=cf_circum, east=cf_circum}}
-- blk0 = FBArray:new{grid=grid, initialState=inflow, label="blk",
-- 		       bcList={west=InFlowBC_Supersonic:new{flowState=inflow},
--                    east=WallBC_NoSlip_FixedT:new{Twall=300, group="loads"},
-- 			       north=OutFlowBC_Simple:new{}}, 
-- 		           nib=1, njb=4}