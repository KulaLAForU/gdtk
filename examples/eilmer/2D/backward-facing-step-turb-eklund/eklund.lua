-- eklund.lua: Turbulent supersonic flow over a backward facing step
-- Dimir Y.X. Pot, Samuel J. Stennett, Wilson Y.K. Chan
-- Ported from Eilmer3, 2018-02-23
-- Reference:
--      Eklund D.R. et al. (1995), Computers & Fluids, 
--      Volume 24, Issue 5, Pages 593-608, 
--      "A comparative computational/experimental investigation 
--      of Mach 2 flow over a rearward-facing step"

config.title = "Eklund's turbulent Mach 2 flow over a backward facing step"
print(config.title)
config.dimensions = 2
config.turbulence_model = "k_omega"
config.viscous = true
config.flux_calculator = 'adaptive'
config.gasdynamic_update_scheme = "classic-rk3"
--
config.max_time = 1.2e-3 -- approx. 5 flow lengths
config.dt_plot =  0.2e-3
config.dt_history = 2.0e-3
config.max_step = 30000000
-- 
config.cfl_value = 0.4	
config.cfl_count = 3
config.stringent_cfl = false -- true is more robust
config.dt_init = 1.0e-14

-- Gas model and initial inflow conditions
nsp, nmodes, gm = setGasModel('ideal-air-gas-model.lua')
p_inf = 35.0e3 	-- Pa
u_inf = 518.0 	-- m/s
T_inf = 167.0	-- K

-- Set up gas state and update thermodynamic transfer coefficients
gas_inf = GasState:new{gm}
gas_inf.p = p_inf; gas_inf.T = T_inf
gm:updateThermoFromPT(gas_inf)
gm:updateSoundSpeed(gas_inf)
gm:updateTransCoeffs(gas_inf)

-- Estimate turbulence quantities for free stream by specifying
-- turbulence intensity and turbulent-to-laminar viscosity ratio
turb_intensity = 0.01
turb_lam_viscosity_ratio = 100.0
tke_inf = 1.5 * (turb_intensity * u_inf)^2
mu_t_inf = turb_lam_viscosity_ratio * gas_inf.mu
omega_inf = gas_inf.rho * tke_inf / mu_t_inf

-- Set up flow state
inflow = FlowState:new{p=p_inf, T=T_inf, velx=u_inf, tke=tke_inf, omega=omega_inf}
print("Inflow Check\n", inflow)

-- Define geometry of the flow domain
x0 = 0.0; y0 = 0.0 -- datum of flow domain (bottom corner of step)
x1 = 0.0345        -- most downstream end of flow domain
y1 = 0.003175      -- height of backward-facing step
y2 = 0.0213        -- max y-height of flow domain

-- To replicate the inflow boundary layers generated by the wind
-- tunnel in the experiments, users can choose to either model the
-- boundary layers with an extended upstream duct length (nozzle),
-- or use a static inflow profile boundary condition.
use_nozzle_flag = false -- default to using static inflow profile BC
x0_noz = -0.08    -- simulated nozzle (extended duct) length
x1_noz = -0.00326 -- end of simulated nozzle length
y1_noz = -x0_noz * math.tan(math.rad(0.5)) + y1
y0_noz = -x1_noz * math.tan(math.rad(0.5)) + y1

-- Build nodes
a0 = Vector3:new{x=x0, y=y0}
b0 = Vector3:new{x=x1, y=y0}
c0 = Vector3:new{x=x1, y=y1}
d0 = Vector3:new{x=x0, y=y1}

a1 = Vector3:new{x=x0, y=y1}
b1 = Vector3:new{x=x1, y=y1}
c1 = Vector3:new{x=x1, y=y2} 
d1 = Vector3:new{x=x0, y=y2}

a2 = Vector3:new{x=x1_noz, y=y0_noz}
b2 = Vector3:new{x=x0, y=y1}
c2 = Vector3:new{x=x0, y=y2} 
d2 = Vector3:new{x=x1_noz, y=y2}

a3 = Vector3:new{x=x0_noz, y=y1_noz}
b3 = Vector3:new{x=x1_noz, y=y0_noz}
c3 = Vector3:new{x=x1_noz, y=y2}
d3 = Vector3:new{x=x0_noz, y=y2}

--  Build blocks
--  ---------------------------------------
--  |            |      |                 | 
--  |  blk_noz   |      |                 |
--  | (optional) | blk2 |      blk1       |
--  |            |      |                 |
--  |            |      |                 |
--  ---------------------------------------
--                      |                 |
--                      |      blk0       |
--                      -------------------

-- Block 0
surf0 = CoonsPatch:new{p00=a0, p10=b0, p11=c0, p01=d0}
cflist0 = { north = RobertsFunction:new{end0=true,end1=false,beta=1.08},
            east  = RobertsFunction:new{end0=true,end1=true, beta=1.01},
            south = RobertsFunction:new{end0=true,end1=false,beta=1.08},
            west  = RobertsFunction:new{end0=true,end1=true, beta=1.01} }
grid0 = StructuredGrid:new{psurface=surf0, niv=141, njv=61, cfList=cflist0}
blk0 = FBArray:new{grid=grid0, nib=5, njb=1, fillCondition=inflow,
                       bcList={north=WallBC_WithSlip:new{}, 
                               east=OutFlowBC_Simple:new{},
                               south=WallBC_NoSlip_Adiabatic:new{}, 
                               west=WallBC_WithSlip:new{}}}

-- Block 1 
surf1 = CoonsPatch:new{p00=a1, p10=b1, p11=c1, p01=d1}
cflist1 = { north = RobertsFunction:new{end0=true,end1=false,beta=1.08},
            east  = RobertsFunction:new{end0=true,end1=true, beta=1.0025},
            south = RobertsFunction:new{end0=true,end1=false,beta=1.08},
            west  = RobertsFunction:new{end0=true,end1=true, beta=1.0025} }
grid1 = StructuredGrid:new{psurface=surf1, niv=141, njv=121, cfList=cflist1}
blk1 = FBArray:new{grid=grid1, nib=5, njb=2, fillCondition=inflow,
                       bcList={north=WallBC_NoSlip_Adiabatic:new{}, 
                               east=OutFlowBC_Simple:new{},
                               south=WallBC_WithSlip:new{}, 
                               west=WallBC_WithSlip:new{}}}

-- Block 2
if use_nozzle_flag then
    west2_BC = WallBC_WithSlip:new{}
else
    west2_BC = InFlowBC_StaticProfile:new{fileName="profile.dat", match="xyA-to-xyA"}
end
--
surf2 = CoonsPatch:new{p00=a2, p10=b2, p11=c2, p01=d2}
cflist2 = { north = RobertsFunction:new{end0=false,end1=true,beta=1.18},
            east  = RobertsFunction:new{end0=true, end1=true,beta=1.0025},
            south = RobertsFunction:new{end0=false,end1=true,beta=1.18}, 
            west  = RobertsFunction:new{end0=true, end1=true,beta=1.0025} }
grid2 = StructuredGrid:new{psurface=surf2, niv=25, njv=121, cfList=cflist2}
blk2 = FBArray:new{grid=grid2, nib=1, njb=2, fillCondition=inflow,
                       bcList={north=WallBC_NoSlip_Adiabatic:new{},
                               east=WallBC_WithSlip:new{},
                               south=WallBC_NoSlip_Adiabatic:new{},
                               west=west2_BC}}
-- Nozzle block
if use_nozzle_flag then
    surf3 = CoonsPatch:new{p00=a3, p10=b3, p11=c3, p01=d3}
    cflist3 = { north = RobertsFunction:new{end0=false,end1=true,beta=1.06},
                east  = RobertsFunction:new{end0=true, end1=true,beta=1.0025},
                south = RobertsFunction:new{end0=false,end1=true,beta=1.06},
                west  = RobertsFunction:new{end0=true, end1=true,beta=1.0025} }
    grid3 = StructuredGrid:new{psurface=surf3, niv=101, njv=121, cfList=cflist3}
    blk_noz = FBArray:new{grid=grid3, nib=1, njb=2, fillCondition=inflow,
                              bcList={north=WallBC_NoSlip_Adiabatic:new{},
                                      east=WallBC_WithSlip:new{},
                                      south=WallBC_NoSlip_Adiabatic:new{},
                                      west=InFlowBC_Supersonic:new{flowCondition=inflow}}}
end

identifyBlockConnections()
