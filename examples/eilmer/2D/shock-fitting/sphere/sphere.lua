-- sphere.lua
-- KD 2016-01-01
-- PJ 2021-08-07 Make same conditions as for cylinder example.
--
config.title = "Sphere in ideal air flow with shock fitting boundary."
print(config.title)
config.dimensions = 2
config.axisymmetric = true
-- problem parameters
u_inf = 2430.0  -- m/s
radius = 1.0    -- m
--
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
initial = FlowState:new{p=100.0e3/3.0, T=200.0, velx=0.0, vely=0.0}
inflow = FlowState:new{p=100.0e3, T=300.0, velx=u_inf, vely=0.0}
--
print "Building grid."
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=-radius, y=0.0}
c = Vector3:new{x=0.0, y=radius}
--
d = Vector3:new{x=-1.5*radius, y=0}
e = Vector3:new{x=-1.5*radius, y=radius}
f = Vector3:new{x=-radius, y=2.0*radius}
g = Vector3:new{x=0.0, y=3.0*radius}
--
psurf = makePatch{north=Line:new{p0=g, p1=c},
		  east=Arc:new{p0=b, p1=c, centre=a},
		  south=Line:new{p0=d, p1=b},
		  west=Bezier:new{points={d, e, f, g}}}
grid = StructuredGrid:new{psurface=psurf, niv=41, njv=41}
--
-- We can leave east and south as slip-walls
blk = FBArray:new{grid=grid, initialState=initial,
                  bcList={west=InFlowBC_ShockFitting:new{flowState=inflow},
                          north=OutFlowBC_Simple:new{}},
                  nib=4, njb=2}
identifyBlockConnections()
--
-- Set a few more config options
body_flow_time = (radius*2)/u_inf
print("body_flow_time=", body_flow_time)
config.flux_calculator = "ausmdv"
config.gasdynamic_update_scheme = "backward_euler"
config.max_time = 20*body_flow_time
config.max_step = 400000
config.cfl_value = 0.5
config.dt_init = 1.0e-7
config.dt_plot = body_flow_time
config.grid_motion = "shock_fitting"
config.shock_fitting_delay = 1*body_flow_time
config.max_invalid_cells = 10
config.adjust_invalid_cell_data = true
config.report_invalid_cells = false
