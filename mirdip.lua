-- set paths to scripts
gkeyllsetup = '/home/space/liang/src/gkeyllscripts/gkeyllsetup.lua'
gkeyllutils = '/home/space/liang/src/gkeyllscripts/gkeyllutils.lua'
-- load helper functions
dofile(gkeyllutils)

----------------
-- PARAMETERS --
----------------
             numMoments = {5, 5}
-- miscellaneous feature switches
hasStairSteppedBoundary = true
     alwaysUseLaxSolver = true
    useIntermediateWave = false
        subTimeStepping = true
              writeDivB = true
              writeDivE = true
                   tEnd = 3600
                nFrames = 60
     Lucee.IsRestarting = false
     Lucee.RestartFrame = -1
-- physical constants
               gasGamma = 5/3
             lightSpeed = Lucee.SpeedOfLight/40
                    mu0 = Lucee.Mu0
               epsilon0 = 1/mu0/(lightSpeed^2)
-- planet radius
                     Re = 6400e3
-- dipole streghth and tilt
                 thetaB = 0
                   phiB = 0
                     B0 = -3.0574e-5
                     Dx = B0*Re^3*math.sin(thetaB)*math.cos(phiB)
                     Dy = B0*Re^3*math.sin(thetaB)*math.sin(phiB)
                     Dz = B0*Re^3*math.cos(thetaB)
                     r1 = 1.5*Re -- zero B field when r < r1
-- radius of inner boundary
                     r0 = 3*Re
-- inflows parameters
                   n_in = 6e6
                  vx_in = 450e3
                  vy_in = 0
                  vz_in = 0
                   p_in = 6e-12
                  Bx_in = 0
                  By_in = 0
                  Bz_in = -5e-9
-- mirdip setup specfications
                r_ramp1 = 12*Re
                r_ramp2 = 14*Re
                stretch = 2
                   xmir = -15*Re
-- kinetic parameters
             mass_ratio = 25
         pressure_ratio = 1
                 d_i_in = 0.5*Re
                ionMass = Lucee.ProtonMass
                elcMass = Lucee.ProtonMass/mass_ratio
                 rho_in = n_in * (elcMass + ionMass)
              ionCharge = ionMass / d_i_in / math.sqrt(mu0*rho_in)
              elcCharge = -ionCharge
                 charge = {elcCharge, ionCharge}
                   mass = {elcMass, ionMass}
-- domain and grid
          numDimensions = 3
           xlo, xup, nx = -20*Re, 20*Re, 40*2
           ylo, yup, ny = -20*Re, 20*Re, 40*2
           zlo, zup, nz = -20*Re, 20*Re, 40*2
                  lower = {xlo, ylo, zlo}
                  upper = {xup, yup, zup}
                  cells = {nx, ny, nz}
-- computational constants
    elcErrorSpeedFactor = 0
    mgnErrorSpeedFactor = 1
     elcErrorDampFactor = 0
     mgnErrorDampFactor = 0
                    cfl = 0.9
-- derived parameters
                  v2_in = vx_in^2 + vy_in^2 + vz_in^2
               rho_e_in = rho_in / (1 + mass_ratio)
                 p_e_in = p_in / (1 + pressure_ratio)
                 u_e_in = p_e_in / (gasGamma-1) + 0.5 * rho_e_in * v2_in
               rho_i_in = rho_in - rho_e_in
                 p_i_in = p_in - p_e_in
                 u_i_in = p_i_in / (gasGamma-1) + 0.5 * rho_i_in * v2_in
                   B_in = math.sqrt(Bx_in^2 + By_in^2 + Bz_in^2)
                pmag_in = B_in^2/2/mu0
                  vA_in = B_in/math.sqrt(mu0*rho_in)
                  cs_in = math.sqrt(gasGamma*p_in/rho_in)
                beta_in = p_in/pmag_in
                 d_i_in = ionMass / ionCharge / math.sqrt(mu0*rho_in)
                  Ex_in = - vy_in*Bz_in + vz_in*By_in
                  Ey_in = - vz_in*Bx_in + vx_in*Bz_in
                  Ez_in = - vx_in*By_in + vy_in*Bx_in
        lightSpeedScale = Lucee.SpeedOfLight / lightSpeed
         elcChargeScale = -Lucee.ElementaryCharge / elcCharge

for _,param in ipairs({
   "gasGamma", "lightSpeed", "mu0",
   "Re", "thetaB", "phiB", "B0", "r0", "r1",
   "n_in", "rho_in", "p_in", "vx_in", "vy_in", "vz_in",
   "Bx_in", "By_in", "Bz_in",
   "vA_in", "cs_in", "beta_in", "d_i_in",
   "mass_ratio", "pressure_ratio", "elcMass", "ionMass", "elcCharge", "ionCharge",
   "lightSpeedScale", "elcChargeScale",
   "xlo", "xup", "nx", "ylo", "yup", "ny", "zlo", "zup", "nz",
}) do
   pretty_log(param)
end

-----------------------
-- INITIAL CONDITION --
-----------------------
function calcRho(x,y,z)
   local rho = rho_in
   return rho
end

function calcP(x,y,z)
   local p = p_in
   return p
end

function calcV(x,y,z)
   local xx = x > 0 and x/stretch or x
   local r = math.sqrt(xx^2 + y^2 + z^2)
   local s = (r-r_ramp1)/(r_ramp2-r_ramp1)
   s = math.max(s, 0)
   s = math.min(s, 1)
   local vx = vx_in * s
   local vy = vy_in * s
   local vz = vz_in * s
   return vx, vy, vz
end

function dipoleB(x, y, z, x0, y0, z0, Dx, Dy, Dz)
   local xx = (x - x0)
   local yy = (y - y0)
   local zz = (z - z0)
   local rr = math.sqrt(xx^2+yy^2+zz^2)
   local Bx = (3*xx*Dx*xx + 3*xx*Dy*yy + 3*xx*Dz*zz - Dx*rr^2) / rr^5
   local By = (3*yy*Dx*xx + 3*yy*Dy*yy + 3*yy*Dz*zz - Dy*rr^2) / rr^5
   local Bz = (3*zz*Dx*xx + 3*zz*Dy*yy + 3*zz*Dz*zz - Dz*rr^2) / rr^5
   if (rr < r1) then
      Bx, By, Bz = 0, 0, 0
   end
   return Bx, By, Bz
end

function staticB(x, y, z)
   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz)
   local Bxs, Bys, Bzs = Bxd0, Byd0, Bzd0
   return Bxs, Bys, Bzs
end

function totalB(x,y,z)
   local Bxd0, Byd0, Bzd0 = dipoleB(x, y, z, 0, 0, 0, Dx, Dy, Dz)
   local Bxdm, Bydm, Bzdm = dipoleB(x, y, z, 2*xmir, 0, 0, -Dx, Dy, Dz)
   local Bxt, Byt, Bzt = Bxd0+Bxdm, Byd0+Bydm, Bzd0+Bzdm
   if (x < xmir) then
      Bxt, Byt, Bzt = 0, 0, 0
   end
   Bxt, Byt, Bzt = Bxt+Bx_in, Byt+By_in, Bzt+Bz_in
   local rr = math.sqrt(x^2+y^2+z^2)
   if (rr < r1) then
      Bxt, Byt, Bzt = 0, 0, 0
   end
   return Bxt, Byt, Bzt
end

function calcB(x,y,z)
   return Bx_in, By_in, Bz_in
end

function init(x,y,z)
   local rho = calcRho(x,y,z)
   local vx,vy,vz = calcV(x,y,z)
   local p = calcP(x,y,z)

   local rho_e = rho / (1 + mass_ratio)
   local rho_i = rho - rho_e
   local rhovx_e = rho_e * vx
   local rhovy_e = rho_e * vy
   local rhovz_e = rho_e * vz
   local rhovx_i = rho_i * vx
   local rhovy_i = rho_i * vy
   local rhovz_i = rho_i * vz
   local p_e = p / (1 + pressure_ratio)
   local p_i = p - p_e

   local u_e = p_e / (gasGamma-1) + 
      0.5 * (rhovx_e^2 + rhovy_e^2 + rhovz_e^2) / rho_e
   local u_i = p_i / (gasGamma-1) + 
      0.5 * (rhovx_i^2 + rhovy_i^2 + rhovz_i^2) / rho_i

   local Bx, By, Bz = calcB(x, y, z)

   local Ex = - vy*Bz + vz*By
   local Ey = - vz*Bx + vx*Bz
   local Ez = - vx*By + vy*Bx

   return rho_e, rhovx_e, rhovy_e, rhovz_e, u_e,
          rho_i, rhovx_i, rhovy_i, rhovz_i, u_i,
          Ex, Ey, Ez, Bx, By, Bz, 0, 0
end

function setStaticField(x, y, z)
   local Exs, Eys, Ezs = 0, 0, 0
   local Bxs, Bys, Bzs = staticB(x, y, z)
   return Exs, Eys, Ezs, Bxs, Bys, Bzs, 0, 0
end

function setInOutField(x, y, z)
   if (math.sqrt(x^2 + y^2 + z^2) < r0) then
      return -1
   else
      return 1
   end
end

------------------------
-- BOUNDARY CONDITION --
------------------------
-- FIXME allow pressure or total energy to float for subsonic inflow BC?
inflowBcElc = BoundaryCondition.Const {
   components = {0, 1, 2, 3, 4},
   values = {rho_e_in, rho_e_in*vx_in, rho_e_in*vy_in, rho_e_in*vz_in, u_e_in},
}
inflowBcIon = BoundaryCondition.Const {
   components = {5, 6, 7, 8, 9},
   values = {rho_i_in, rho_i_in*vx_in, rho_i_in*vy_in, rho_i_in*vz_in, u_i_in},
}
inflowBcE = BoundaryCondition.Const {
   components = {10, 11, 12},
   values = { Ex_in, Ey_in, Ez_in },
}
bcInflowB = BoundaryCondition.Const {
   components = {13, 14, 15},
   values = { Bx_in, By_in, Bz_in },
}
inflowBcPot = BoundaryCondition.Copy { components = {16, 17} }

function createInflowBc()
   return Updater.Bc3D {
      onGrid = grid,
      boundaryConditions = { 
         inflowBcElc, inflowBcIon, inflowBcE, inflowBcB, inflowBcPot
      },
      dir = 0,
      edge = 'lower',
   }
end

bcElcCopy = BoundaryCondition.Copy { components = {0, 4} }
bcElcWall = BoundaryCondition.ZeroNormal { components = {1, 2, 3} }
bcIonCopy = BoundaryCondition.Copy { components = {5, 9} }
bcIonWall = BoundaryCondition.ZeroNormal { components = {6, 7, 8} }
bcEConduct = BoundaryCondition.ZeroTangent { components = {10, 11, 12} }
bcBConduct = BoundaryCondition.ZeroNormal { components = {13, 14, 15} }
bcPotCopy = BoundaryCondition.Copy { components = {16, 17}, fact = {1, 1} }

function createInnerBc()
   return Updater.StairSteppedBc3D {
      onGrid = grid,
      inOutField = inOutField,
      boundaryConditions = {
         bcElcCopy, bcElcWall, bcIonCopy, bcIonWall,
         bcEConduct, bcBConduct, bcPotCopy,
      }
   }
end

function applyBc(myQ, tCurr, myDt, dir)
   if not inflowBc then
      inflowBc = createInflowBc()
   end
   runUpdater(inflowBc, tCurr, myDt, nil, {myQ}, nil)
   myQ:applyCopyBc(0, "upper")
   myQ:applyCopyBc(1, "lower")
   myQ:applyCopyBc(1, "upper")
   myQ:applyCopyBc(2, "lower")
   myQ:applyCopyBc(2, "upper")
   if not innerBc then
      innerBc = createInnerBc()
   end
   runUpdater(innerBc, tCurr, myDt, nil, {myQ}, dir)
   myQ:sync()
end

-----------
-- SETUP --
-----------
dofile(gkeyllsetup)

-------------------
-- TIME-STEPPING --
-------------------
tStart = 0
initDt = 100
runSimulation(tStart, tEnd, nFrames, initDt)
