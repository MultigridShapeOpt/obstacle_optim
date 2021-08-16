  
  -- * Author: Jose Pinzon
   -- Source: https://github.com/MultigridShapeOpt
  -- *
  -- * This file is a part of the FluidOptim UG4 plugin under development at 
  -- * the Research Group Approximation and Optimization, Hamburg University
  -- * and as part of the project SENSUS (LFF-GK11).
  -- *
  -- * This library is free software; you can redistribute it and/or
  -- * modify it under the terms of the GNU Lesser General Public
  -- * License as published by the Free Software Foundation; either
  -- * version 2.1 of the License, or (at your option) any later version.
  -- *
  -- * This library is distributed in the hope that it will be useful,
  -- * but WITHOUT ANY WARRANTY; without even the implied warranty of
  -- * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  -- * Lesser General Public License for more details.


PluginRequired("FluidOptim")
-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")
ug_load_script("util/conv_rates_static.lua")
ug_load_script("util/solver_util.lua")
ug_load_script("util/load_balancing_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("obstacle_optim_3d_util.lua")
ug_load_script("ns_util.lua")

ProfileLUA(true)
PrintBuildConfiguration()

package.path = package.path .. ";../../apps/obstacle_optim/lbfgs/?.lua"--Jose line
--package.path = package.path .. ";../lbfgs/?.lua"--Prof Siebenborn line
local deque = require("deque")
local vbfgs = require("vbfgs")
local lbfgs = require("lbfgs")
--Constants
dim = 3
numPreRefs=0
vorder=1
porder=1
d_initial_value=1.0;
deforder=1
diameter=6


--Parameters
numRefs = util.GetParamNumber("-numRefs", 1)--refinements after grid partitioning/distribution
visc = util.GetParamNumber("-visc",1.0e0)--medium viscosity
stab = util.GetParamNumber("-stab", 0.0)--stabilization for navier stokes
threshold = util.GetParamNumber("-threshold", 1.0e-1)--determinant threshold
beta = util.GetParamNumber("-beta", 100.0e0)--beta value, must be high
alpha =  util.GetParamNumber("-alpha", 5.0e-3)--regularization parameter alpha, must be small
nu_geom = util.GetParamNumber("-nu_geom", 60.0e0)--geometrical penalty, must be high
numSteps= util.GetParamNumber("-numSteps", 5000)
stepSizeBnd= util.GetParamNumber("-stepSizeBnd", 1e-1)
stepSizeExt= util.GetParamNumber("-stepSizeExt", 1e-1)
updateParam= util.GetParamNumber("-updateParam", 5)
ext_upper = util.GetParamNumber("-ext_upper",60.00)
ext_lower = util.GetParamNumber("-ext_lower", 0.0)
nu_ext = util.GetParamNumber("-nu_ext", (ext_upper-ext_lower)/2)--extension factor
chi = util.GetParamNumber("-chi", 1.0e-8)--penalization on the mean of the extension factor, increase for smoothness
limit_memory = util.GetParamNumber("-limit_memory", 3)
control_tol = util.GetParamNumber("-control_tol", 1e0)
geom_penalty_inc = util.GetParamNumber("-geom_penalty_inc", 1.5e-0)
geom_constraints_tol = util.GetParamNumber("-geom_constraints_tol", 1e-2)
multiplier_scaling = util.GetParamNumber("-multiplier_scaling", 5e-1)
stabType = util.GetParamNumber("-stabType", 80.0)
ext_overwrite=util.GetParamNumber("-ext_overwrite",ext_upper)
printInterval=util.GetParamNumber("-printInterval",1)--1 print all
vorder=util.GetParamNumber("-vorder",1)--1 print all

bOutput  = util.GetParamBool("-bOutput",false)--output VTK visualization files, all discretizations
bOutputPlot = util.GetParamBool("-bOutputPlot",false)--plots with gnuplot
bDoNothing = util.GetParamBool("-bDoNothing",true)--do nothing on flow outlet
bDeformationVolumeOutput = util.GetParamBool("-bDeformationVolumeOutput",true)--output deformation vector field
bDeformationSurfaceOutput = util.GetParamBool("-bDeformationSurfaceOutput",true)--output deformation vector field
bNavierStokesOutput = util.GetParamBool("-bNavierStokesOutput",false)--output deformation vector field
bAdjointFlowsOutput = util.GetParamBool("-bAdjointFlowsOutput",false)--output deformation vector field
bAdjointDeformationsOutput = util.GetParamBool("-bAdjointDeformationsOutput",false)--output deformation vector field
bExtensionEqOutput = util.GetParamBool("-bExtensionEqOutput",true)--output K scalar field
bDesignEqOutput = util.GetParamBool("-bDesignEqOutput",true)--output G scalar field
bLoadPreviousSteps = util.GetParamBool("-bLoadPreviousSteps",false)
bPrintSingleStep = util.GetParamBool("-bPrintSingleStep", false)--override vtk output step
gridName = util.GetParam("-grid", "./grids/ns.ugx")

--Variable names, still have to manually change boundary condition definitions

if dim == 3 then 
deformationNames="u1,u2,u3"
flowNames="v1,v2,v3"
pressureName="p"
adjointFlowNames="q1,q2,q3"
adjointPressure="h"
adjointDeformationNames="l1,l2,l3"
designEquationNames = "g"
ExtensionEquationNames = "k"
elseif dim ==2 then 
deformationNames="u1,u2"
flowNames="v1,v2"
pressureName="p"
adjointFlowNames="q1,q2"
adjointPressure="h"
adjointDeformationNames="l1,l2"
designEquationNames = "g"
ExtensionEquationNames = "k"
end
--Output file names
deformationOutputFile="vtk_deformation"
flowOutputFile="vtk_flows"
adjointFlowOutputFile="vtk_adjointFlows"
adjointDeformationOutputFile="vtk_adjointDeformation"
designEquationOutputFile="vtk_designEquation"
ExtensionEquationOutputFile="vtk_extensionfactor"
--Vector names
deformationField="deformation"
nodalFlows="flows"
adjointNodalFlows="adjointFlows"
adjointDeformationField="adjoint_deformation"
designEquation="control_variable"
ExtensionFactor="extension_factor"
--Remarks
print("USAGE REMARKS: ")
print("If the function TransformDomainByDisplacement is used (deformation of domain), the scheme collapses")
print("GRID FUNCTIONS ARE:")
print("w: deformation field")
print("v: nodal flow vector, pressure")
print("q: adjoint nodal flow vector, pressure")
print("l: adjoint deformations")
print("d: boundary control variable field")
print("g: gradient of boundary control")
print("k: gradient of extension factor")

print("THE PARAMETERS USED FOR EXECUTION ARE: ")
print("grid: "..gridName)
print("velocity order "..vorder)
print("pressure order "..porder)
print("numPreRefs:   ".. numPreRefs)
print("numRefs:      ".. numRefs)
print("numSteps:     ".. numSteps)
print("lagrange. update: ".. updateParam)
print("stabilization:".. stab)
print("viscosity:    ".. visc)
print("beta:         ".. beta)
print("alpha:        ".. alpha)
print("nu_geom:     " .. nu_geom)
print("chi:".. chi)
print("ext_lower:    ".. ext_lower)
print("ext_upper:         ".. ext_upper)
print("initial extension factor:     "..nu_ext)
print("determinant threshold:    ".. threshold)
print("step size on boundary control:	"..stepSizeBnd)
print("step size on extension factor:	"..stepSizeExt)
print("bfgs memory:   "..limit_memory)
print("geom_penalty_inc:   "..geom_penalty_inc)
print("control_tol:   "..control_tol)
print("geom_constraints_tol:   "..geom_constraints_tol)
print("stabilization technique "..stabType)

-- initialize ug with the world dimension and the algebra type
InitUG(dim, AlgebraType("CPU", 1));

-- load grid into domain
dom = Domain()
LoadDomain(dom, gridName)
--dom = util.CreateDomain(gridName, 0, {})
print("Loaded domain from " .. gridName)

--[[
local refiner =  GlobalDomainRefiner(dom)
for i=1,numPreRefs do refiner:refine(); end
if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, false) == false then
	--print("Error while Distributing Grid. Aborting."); exit();
end
for i=numPreRefs+1,numRefs do refiner:refine(); end
--]]

-- Refine the domain (redistribution is handled internally for parallel runs)
--print("refining...")

----[[
-- This balancing setup makes sense for structured grids with uniform refinement
balancerDesc = {
--[[
                partitioner = {
                        name = "dynamicBisection",
                        verbose = false,
                        enableXCuts = true,
                        enableYCuts = true,
                        enableZCuts = true,
                        longestSplitAxis = false,
                        clusteredSiblings = true,
                        balanceThreshold = 0.9,
                        numSplitImprovements = 10
                },
--]]
----[[
                partitioner = {
                        name = "staticBisection",
                        verbose = false,
                        enableXCuts = true,
                        enableYCuts = true,
                        enableZCuts = true,
                        longestSplitAxis = false,
                        clusteredSiblings = true,
                        balanceThreshold = 0.9,
                        numSplitImprovements = 10
                },
--]]

		hierarchy = {
			name 						= "noRedists",
			minElemsPerProcPerLevel		= 4,--better for coarser grids
			maxRedistProcs				= 120,
			
			{-- levels 0 to 0
			upperLvl = 0,
			maxProcs = 1
            },
		},
	}

util.refinement.CreateRegularHierarchy(dom, numRefs, false, balancerDesc)
--]]
print(dom:domain_info():to_string())

print("DEFORMATION STATE EQUATION: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
deformationState_ApproxSpace = ApproximationSpace(dom)
--deformationState_ApproxSpace:add_fct(deformationNames, "Lagrange", 1)
deformationState_ApproxSpace:add_fct("u1", "Lagrange", 1)
deformationState_ApproxSpace:add_fct("u2", "Lagrange", 1)
if dim == 3 then
	deformationState_ApproxSpace:add_fct("u3", "Lagrange", 1)
end
deformationState_ApproxSpace:init_levels()
deformationState_ApproxSpace:init_top_surface()
print("DEFORMATION STATE EQUATION: Approx. Space:")
deformationState_ApproxSpace:print_statistic()

deformationState_elemDisc=DeformationStateEquation(deformationNames,"outer")
deformationState_elemDisc:set_symmetric(true)
--nonlinear settings
deformationState_elemDisc:set_nonlinear(true)
deformationState_elemDisc:set_picard(false)
deformationState_elemDisc:set_extension_factor(nu_ext)

boundaryControl_elemDisc=SurfaceBoundaryControl(deformationNames,"obstacle_surface")
boundaryControl_elemDisc:set_boundary_control(0.0)


deformationState_Dirich = DirichletBoundary()
--constraints on x
deformationState_Dirich:add(0.0,"u1","wall")
deformationState_Dirich:add(0.0,"u1","inlet")
deformationState_Dirich:add(0.0,"u1","outlet")
--constraints on y
deformationState_Dirich:add(0.0,"u2","wall")
deformationState_Dirich:add(0.0,"u2","inlet")
deformationState_Dirich:add(0.0,"u2","outlet")
if dim == 3 then 
--constraints on z
deformationState_Dirich:add(0.0,"u3","wall")
deformationState_Dirich:add(0.0,"u3","inlet")
deformationState_Dirich:add(0.0,"u3","outlet")
end
--Domain Discretization

deformationState_DomainDisc = DomainDiscretization(deformationState_ApproxSpace)
deformationState_DomainDisc:add(deformationState_elemDisc)
deformationState_DomainDisc:add(deformationState_Dirich)
deformationState_DomainDisc:add(boundaryControl_elemDisc)
print("DEFORMATION STATE EQUATION: create GridFunctions, Matrix Operator, and GlobalGridFunctionGradientDatas")
w = GridFunction(deformationState_ApproxSpace);w:set(0.0)
deformationState_DomainDisc:adjust_solution(w)
w1_gradient_global=GlobalGridFunctionGradientData(w,"u1")
w2_gradient_global=GlobalGridFunctionGradientData(w,"u2")
w1_value_global=GlobalGridFunctionNumberData(w,"u1")
w2_value_global=GlobalGridFunctionNumberData(w,"u2")

if dim == 3 then
w3_gradient_global=GlobalGridFunctionGradientData(w,"u3")
w3_value_global=GlobalGridFunctionNumberData(w,"u3")
end

deformation_Op = AssembledOperator()
deformation_Op:set_discretization(deformationState_DomainDisc)
print("DEFORMATION STATE EQUATION: all set")

--NAVIER STOKES 
function CreateNSSolver(approxSpace, p, linsolver, discretization)

	local base = LU()

	local smoother = nil
	if porder == vorder then
		smoother = ILU()
		smoother:set_damp(1.0)
	else
		--smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
		smoother = ILU()
		smoother:set_damp(0.7)
	end

	local gmg = util.gmg.create(approxSpace, smoother, 3, 3,
							 "V", base, 0, false, discretization)
	gmg:set_rap(true)
	gmg:set_gathered_base_solver_if_ambiguous(false)
	transfer = StdTransfer()
	gmg:set_transfer(transfer)
	
	--gmg:add_prolongation_post_process(AverageComponent("p"))
	-- transfer:enable_p1_lagrange_optimization(false)
	
	local solver = util.solver.create(linsolver, gmg)	
	solver:set_convergence_check(ConvCheck(10000, 1e-12, 1e-2, true))	
	
	local convCheck = ConvCheck(500, 1e-8, 1e-99, true)	
	local newtonSolver = NewtonSolver()
	newtonSolver:set_linear_solver(solver)
	newtonSolver:set_convergence_check(convCheck)
	lineSearch=StandardLineSearch(50, 1.0, 0.9, true, false)
	lineSearch:set_verbose(false)
	--newtonSolver:set_line_search(StandardLineSearch(50, 1.0, 0.9, true, false))
	newtonSolver:set_line_search(lineSearch)
	--newtonSolver:line_search():set_verbose(false)
	--newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
	
	return newtonSolver
end

function CreateAdjointFlowSolver(approxSpace, p, linsolver, discretization)

	local base = LU()
	
	local smoother = nil
	if porder == vorder then
		smoother = ILU()
		smoother:set_damp(0.7)
	else
		--smoother = ComponentGaussSeidel(0.1, {"p"}, {0}, {1})
		smoother = ILU()
		smoother:set_damp(1.0)
	end
	
	local gmg = util.gmg.create(approxSpace, smoother, 3, 3,
							 "V", base, 0, false, discretization)
	gmg:set_rap(true)
	gmg:set_gathered_base_solver_if_ambiguous(false)
	transfer = StdTransfer()
	gmg:set_transfer(transfer)
	
	--gmg:add_prolongation_post_process(AverageComponent("p"))
	-- transfer:enable_p1_lagrange_optimization(false)

	local solver = util.solver.create(linsolver, gmg)	
	solver:set_convergence_check(ConvCheck(10000, 1e-10, 1e-10, true))	

	
	return solver
end

--DIRICHLET BOUNDARY VALUES AT INLET
function SquareInletVelocities(x,y,z,t)
	local s1=math.abs(y)*math.pi/diameter
	local s2=math.abs(z)*math.pi/diameter
	return math.max(0.0,math.cos(s1)*math.cos(s2))
	--return 0
end

function InletVelocities(x,y,z,t)
	local s=math.sqrt(y*y+z*z)*math.pi/diameter
	return math.cos(s)
	--return 0
end

print("NAVIER STOKES: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
NavierStokes_ApproxSpace = ApproximationSpace(dom)
NavierStokes_ApproxSpace:add_fct(flowNames,"Lagrange", vorder)
NavierStokes_ApproxSpace:add_fct(pressureName,"Lagrange",porder) 
NavierStokes_ApproxSpace:init_levels()
NavierStokes_ApproxSpace:init_top_surface()

print("NAVIER STOKES: Approx. Space:")
NavierStokes_ApproxSpace:print_statistic()
NavierStokes_ElemDisc = nil
if dim == 2 then 
NavierStokes_ElemDisc = TransformedNavierStokes("v1,v2,p ", "outer")
end
if dim == 3 then 
NavierStokes_ElemDisc = TransformedNavierStokes("v1,v2,v3,p ", "outer")
end
--useful settings
NavierStokes_ElemDisc:set_nonlinear(true)
NavierStokes_ElemDisc:set_picard(false)
NavierStokes_ElemDisc:no_deformation(false)

NavierStokes_ElemDisc:set_kinematic_viscosity(visc)
NavierStokes_ElemDisc:set_stabilization(stab)
NavierStokes_ElemDisc:set_stabilization_type(stabType)

NavierStokes_ElemDisc:set_deformation_vector_d1(w1_gradient_global)
NavierStokes_ElemDisc:set_deformation_vector_d2(w2_gradient_global)
if dim == 3 then 
NavierStokes_ElemDisc:set_deformation_vector_d3(w3_gradient_global)
end
--************BOUNDARY CONDITIONS*********--
NavierStokes_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
NavierStokes_Dirich:add("SquareInletVelocities","v1","inlet")
--tns_dirichletBND:add(1,"u","inlet")
NavierStokes_Dirich:add(0,"v2","inlet")
--************WALL BOUNDARY**************--
NavierStokes_Dirich:add(0,"v1","wall")
NavierStokes_Dirich:add(0,"v2","wall")
--************OBSTACLE SURFACE BOUNDARY**************--
NavierStokes_Dirich:add(0,"v1","obstacle_surface")
NavierStokes_Dirich:add(0,"v2","obstacle_surface")
if dim == 3 then 
--************INLET BOUNDARY**************--
NavierStokes_Dirich:add(0,"v3","inlet")
--************WALL BOUNDARY**************--
NavierStokes_Dirich:add(0,"v3","wall")
--************OBSTACLE SURFACE BOUNDARY**************--
NavierStokes_Dirich:add(0,"v3","obstacle_surface")
end

if not bDoNothing then
	print("Output flows set MIND YOU!!!!!!!!!!!!!!!!!!!!!")
	NavierStokes_Dirich:add("InletVelocities","v1","outlet")
	NavierStokes_Dirich:add(0,"v2","outlet")
	NavierStokes_Dirich:add(0,"v3","outlet")
end
-- Domain Discretization 
NavierStokes_DomainDisc = DomainDiscretization(NavierStokes_ApproxSpace)
NavierStokes_DomainDisc:add(NavierStokes_ElemDisc)
NavierStokes_DomainDisc:add(NavierStokes_Dirich)
print("NAVIER STOKES: create GridFunctions, Matrix Operator, and GlobalGridFunctionGradientDatas")
v = GridFunction(NavierStokes_ApproxSpace);v:set(0.0)
NavierStokes_DomainDisc:adjust_solution(v)
v1_gradient_global=GlobalGridFunctionGradientData(v,"v1")
v2_gradient_global=GlobalGridFunctionGradientData(v,"v2")
v1_value_global=GlobalGridFunctionNumberData(v,"v1")
v2_value_global=GlobalGridFunctionNumberData(v,"v2")
if dim == 3 then 
v3_gradient_global=GlobalGridFunctionGradientData(v,"v3")
v3_value_global=GlobalGridFunctionNumberData(v,"v3")
end
---Pressure
p_value_global=GlobalGridFunctionNumberData(v,"p")
--Nonlinear Matrix Operator
navier_Op = AssembledOperator()
navier_Op:set_discretization(NavierStokes_DomainDisc)
print("NAVIER STOKES: all set")


print("ADJOINT FLOW: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
AdjointFlow_ApproxSpace = ApproximationSpace(dom)
AdjointFlow_ApproxSpace:add_fct(adjointFlowNames,"Lagrange",vorder)
AdjointFlow_ApproxSpace:add_fct(adjointPressure,"Lagrange",porder) 
AdjointFlow_ApproxSpace:init_levels()
AdjointFlow_ApproxSpace:init_top_surface()
print("ADJOINT FLOW : Approx. Space:")
AdjointFlow_ApproxSpace:print_statistic()

AdjointFlow_ElemDisc = nil
if dim == 2 then 
AdjointFlow_ElemDisc = AdjointSystem("q1,q2, h ", "outer")
end
if dim == 3 then
AdjointFlow_ElemDisc = AdjointSystem("q1,q2,q3, h ", "outer")
end
AdjointFlow_ElemDisc:set_kinematic_viscosity(visc)
AdjointFlow_ElemDisc:set_stabilization(stab)
AdjointFlow_ElemDisc:set_stabilization_type(stabType)

--Set Imports
--DEFORMATION GRADIENT
AdjointFlow_ElemDisc:set_deformation_vector_d1(w1_gradient_global)
AdjointFlow_ElemDisc:set_deformation_vector_d2(w2_gradient_global)

--VELOCITY GRADIENT
AdjointFlow_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
AdjointFlow_ElemDisc:set_velocity_vector_d2(v2_gradient_global);

--VELOCITY VECTOR
AdjointFlow_ElemDisc:set_velocity_d1(v1_value_global)
AdjointFlow_ElemDisc:set_velocity_d2(v2_value_global)

if dim == 3 then 
--DEFORMATION GRADIENT
AdjointFlow_ElemDisc:set_deformation_vector_d3(w3_gradient_global)
--VELOCITY GRADIENT
AdjointFlow_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
--VELOCITY VECTOR
AdjointFlow_ElemDisc:set_velocity_d3(v3_value_global)
end
--************BOUNDARY CONDITIONS*********--
AdjointFlow_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","inlet")
AdjointFlow_Dirich:add(0,"q2","inlet")
--************WALL BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","wall")
AdjointFlow_Dirich:add(0,"q2","wall")
--************OBSTACLE SURFACE BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q1","obstacle_surface")
AdjointFlow_Dirich:add(0,"q2","obstacle_surface")
--TODO:see if do nothing is ok, might be adjustable think about it...
--************OUTLET BOUNDARY**************--
--AdjointFlow_Dirich:add(0,"q1","outlet")
--AdjointFlow_Dirich:add(0,"q2","outlet")

if dim == 3 then 
--************INLET BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q3","inlet")
--************WALL BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q3","wall")
--************OBSTACLE SURFACE BOUNDARY**************--
AdjointFlow_Dirich:add(0,"q3","obstacle_surface")
end

-- Domain Discretization 
AdjointFlow_DomainDisc = DomainDiscretization(AdjointFlow_ApproxSpace)
AdjointFlow_DomainDisc:add(AdjointFlow_ElemDisc)
AdjointFlow_DomainDisc:add(AdjointFlow_Dirich)

print("ADJOINT FLOW: create GridFunctions and GlobalGridFunctionGradientDatas")
q = GridFunction(AdjointFlow_ApproxSpace);q:set(0.0)
r_q = GridFunction(AdjointFlow_ApproxSpace);r_q:set(0.0)
AdjointFlow_DomainDisc:adjust_solution(q)
q1_gradient_global=GlobalGridFunctionGradientData(q,"q1")
q2_gradient_global=GlobalGridFunctionGradientData(q,"q2")
q1_value_global=GlobalGridFunctionNumberData(q,"q1")
q2_value_global=GlobalGridFunctionNumberData(q,"q2")
---Pressure
h_value_global=GlobalGridFunctionNumberData(q,"h")
A = AssembledLinearOperator(AdjointFlow_DomainDisc)
print("ADJOINT FLOW SYSTEM: all set:")

if dim == 3 then 
q3_gradient_global=GlobalGridFunctionGradientData(q,"q3")
q3_value_global=GlobalGridFunctionNumberData(q,"q3")
end

ReferenceVolume = VolumeDefect(w,0, "outer", deformationNames,4,false,1,false)--used as a constant
vBarycenterDefect = {0,0,0}
volumeDefect = 0
LambdaVolume = 0
LambdaBarycenter={0,0,0}

print("ADJOINT DEFORMATION: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
AdjointDeformation_ApproxSpace = ApproximationSpace(dom)
AdjointDeformation_ApproxSpace:add_fct(adjointDeformationNames,"Lagrange",1)
AdjointDeformation_ApproxSpace:init_levels()
AdjointDeformation_ApproxSpace:init_top_surface()
print("ADJOINT DEFORMATION : Approx. Space:")
AdjointDeformation_ApproxSpace:print_statistic()

AdjointDeformation_ElemDisc = DeformationAdjointSystem(adjointDeformationNames,"outer") 
AdjointDeformation_ElemDisc:set_kinematic_viscosity(visc)
AdjointDeformation_ElemDisc:set_extension_factor(nu_ext)
AdjointDeformation_ElemDisc:set_determinant_threshold(threshold)
AdjointDeformation_ElemDisc:set_beta(beta)
AdjointDeformation_ElemDisc:set_nu_penalty(nu_geom)
AdjointDeformation_ElemDisc:set_stabilization_scale(stabType)
--Defects and Lagrange Multipliers Initial Set
AdjointDeformation_ElemDisc:set_barycenter_defect(vBarycenterDefect[1],vBarycenterDefect[2],vBarycenterDefect[3])
AdjointDeformation_ElemDisc:set_volume_defect(volumeDefect)
AdjointDeformation_ElemDisc:set_lambda_vol(LambdaVolume)
AdjointDeformation_ElemDisc:set_lambda_barycenter(LambdaBarycenter[1],LambdaBarycenter[2],LambdaBarycenter[3])

--Set Imports
--DEFORMATION VECTOR
AdjointDeformation_ElemDisc:set_deformation_d1(w1_value_global)
AdjointDeformation_ElemDisc:set_deformation_d2(w2_value_global)
--DEFORMATION GRADIENT
AdjointDeformation_ElemDisc:set_deformation_vector_d1(w1_gradient_global)
AdjointDeformation_ElemDisc:set_deformation_vector_d2(w2_gradient_global)

--VELOCITY GRADIENT
AdjointDeformation_ElemDisc:set_velocity_vector_d1(v1_gradient_global);
AdjointDeformation_ElemDisc:set_velocity_vector_d2(v2_gradient_global);
--VELOCITY VECTOR
AdjointDeformation_ElemDisc:set_velocity_d1(v1_value_global)
AdjointDeformation_ElemDisc:set_velocity_d2(v2_value_global)
--PRESSURE
AdjointDeformation_ElemDisc:set_pressure(p_value_global)

--ADJOINT VELOCITY VECTOR
AdjointDeformation_ElemDisc:set_adjoint_velocity_d1(q1_value_global)
AdjointDeformation_ElemDisc:set_adjoint_velocity_d2(q2_value_global)
--ADJOINT VELOCITY GRADIENT
AdjointDeformation_ElemDisc:set_adjoint_velocity_vector_d1(q1_gradient_global)
AdjointDeformation_ElemDisc:set_adjoint_velocity_vector_d2(q2_gradient_global)
--ADJOINT PRESSURE
AdjointDeformation_ElemDisc:set_adjoint_pressure(h_value_global)

if dim == 3 then 
--DEFORMATION VECTOR
AdjointDeformation_ElemDisc:set_deformation_d3(w3_value_global)
--DEFORMATION GRADIENT
AdjointDeformation_ElemDisc:set_deformation_vector_d3(w3_gradient_global)
--VELOCITY GRADIENT
AdjointDeformation_ElemDisc:set_velocity_vector_d3(v3_gradient_global);
--VELOCITY VECTOR
AdjointDeformation_ElemDisc:set_velocity_d3(v3_value_global)
--ADJOINT VELOCITY VECTOR
AdjointDeformation_ElemDisc:set_adjoint_velocity_d3(q3_value_global)
--ADJOINT VELOCITY GRADIENT
AdjointDeformation_ElemDisc:set_adjoint_velocity_vector_d3(q3_gradient_global)
end

--************BOUNDARY CONDITIONS*********--
AdjointDeformation_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
AdjointDeformation_Dirich:add(0,"l1","inlet")
AdjointDeformation_Dirich:add(0,"l2","inlet")
--************WALL BOUNDARY**************--
AdjointDeformation_Dirich:add(0,"l1","wall")
AdjointDeformation_Dirich:add(0,"l2","wall")
--************OUTLET BOUNDARY**************--
AdjointDeformation_Dirich:add(0,"l1","outlet")
AdjointDeformation_Dirich:add(0,"l2","outlet")

if dim == 3 then 
--************INLET BOUNDARY**************--
AdjointDeformation_Dirich:add(0,"l3","inlet")
--************WALL BOUNDARY**************--
AdjointDeformation_Dirich:add(0,"l3","wall")
--************OUTLET BOUNDARY**************--
AdjointDeformation_Dirich:add(0,"l3","outlet")
end

-- Domain Discretization 
AdjointDeformation_DomainDisc = DomainDiscretization(AdjointDeformation_ApproxSpace)
AdjointDeformation_DomainDisc:add(AdjointDeformation_ElemDisc)
AdjointDeformation_DomainDisc:add(AdjointDeformation_Dirich)

print("ADJOINT DEFORMATION: create GridFunctions and GlobalGridFunctionGradientDatas")
l = GridFunction(AdjointDeformation_ApproxSpace);l:set(0.0)
r_l = GridFunction(AdjointDeformation_ApproxSpace);r_l:set(0.0)
AdjointDeformation_DomainDisc:adjust_solution(l)
l1_value_global=GlobalGridFunctionNumberData(l,"l1")
l2_value_global=GlobalGridFunctionNumberData(l,"l2")
if dim == 3 then 
l3_value_global=GlobalGridFunctionNumberData(l,"l3")
end
A_l = AssembledLinearOperator(AdjointDeformation_DomainDisc)
print("ADJOINT DEFORMATION: all set:")


print("SURFACE EQUATION: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
SurfaceDesignEquation_ApproxSpace = ApproximationSpace(dom)
SurfaceDesignEquation_ApproxSpace:add_fct("g","Lagrange",1)
SurfaceDesignEquation_ApproxSpace:init_levels()
SurfaceDesignEquation_ApproxSpace:init_top_surface()
print("Surface DESIGN EQUATION : Approx. Space:")
SurfaceDesignEquation_ApproxSpace:print_statistic()

--Elem/DomainDiscretization
SurfaceDesignEquation_ElemDisc = SurfaceDesignEquation("g","obstacle_surface")
SurfaceRHS_ElemDisc = SurfaceRHS("g","obstacle_surface")
SurfaceRHS_ElemDisc:set_boundary_control(0)
SurfaceRHS_ElemDisc:set_adjoint_deformation_d1(l1_value_global)
SurfaceRHS_ElemDisc:set_adjoint_deformation_d2(l2_value_global)
if dim == 3 then 
SurfaceRHS_ElemDisc:set_adjoint_deformation_d3(l3_value_global)
end
SurfaceRHS_ElemDisc:set_alpha(alpha)

--************BOUNDARY CONDITIONS*********--
SurfaceDesignEquation_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
SurfaceDesignEquation_Dirich:add(0,"g","inlet")
--************WALL BOUNDARY**************--
SurfaceDesignEquation_Dirich:add(0,"g","wall")
--************OUTLET BOUNDARY**************--
SurfaceDesignEquation_Dirich:add(0,"g","outlet")
SurfaceDesignEquation_Dirich:add(0,"g","outer")

-- Domain Discretization 
SurfaceDesignEquation_DomainDisc = DomainDiscretization(SurfaceDesignEquation_ApproxSpace)
SurfaceDesignEquation_DomainDisc:add(SurfaceDesignEquation_ElemDisc)
SurfaceDesignEquation_DomainDisc:add(SurfaceRHS_ElemDisc)
SurfaceDesignEquation_DomainDisc:add(SurfaceDesignEquation_Dirich)
print("SURFACE DESIGN EQUATION: create GridFunctions and GlobalGridFunctionGradientDatas")
g = GridFunction(SurfaceDesignEquation_ApproxSpace);g:set(0.0)
d=GridFunction(SurfaceDesignEquation_ApproxSpace);d:set(0.0)
r_g = GridFunction(SurfaceDesignEquation_ApproxSpace);r_g:set(0.0)
SurfaceDesignEquation_DomainDisc:adjust_solution(g)
A_g = AssembledLinearOperator(SurfaceDesignEquation_DomainDisc)
d_global=GlobalGridFunctionNumberData(d,"g")
print("SURFACE DESIGN EQUATION: all set:")

print("EXTENSION EQUATION: initialize approxSpace, Elem/DomainDiscretization, boundary conds.")
-- set up approximation space
ExtensionEquation_ApproxSpace = ApproximationSpace(dom)
ExtensionEquation_ApproxSpace:add_fct("k","Lagrange",1)
ExtensionEquation_ApproxSpace:init_levels()
ExtensionEquation_ApproxSpace:init_top_surface()
print("EXTENSION EQUATION : Approx. Space:")
ExtensionEquation_ApproxSpace:print_statistic()

--Elem/DomainDiscretization
ExtensionEquation_ElemDisc = ExtensionEquation("k","outer")
ExtensionEquation_ElemDisc:set_extension_factor(nu_ext)
ExtensionEquation_ElemDisc:set_extension_lower(ext_lower)
ExtensionEquation_ElemDisc:set_extension_upper(ext_upper)
ExtensionEquation_ElemDisc:set_chi(chi)
--ExtensionEquation_ElemDisc:set_quad_order(1)
--Set Imports
--DEFORMATION VECTOR
ExtensionEquation_ElemDisc:set_deformation_d1(w1_value_global)
ExtensionEquation_ElemDisc:set_deformation_d2(w2_value_global)
--DEFORMATION GRADIENT
ExtensionEquation_ElemDisc:set_deformation_vector_d1(w1_gradient_global)
ExtensionEquation_ElemDisc:set_deformation_vector_d2(w2_gradient_global)
--ADJOINT DEFORMATION
ExtensionEquation_ElemDisc:set_adjoint_deformation_d1(l1_value_global)
ExtensionEquation_ElemDisc:set_adjoint_deformation_d2(l2_value_global)

if dim == 3 then 
--DEFORMATION VECTOR
ExtensionEquation_ElemDisc:set_deformation_d3(w3_value_global)
--DEFORMATION GRADIENT
ExtensionEquation_ElemDisc:set_deformation_vector_d3(w3_gradient_global)
--ADJOINT DEFORMATION
ExtensionEquation_ElemDisc:set_adjoint_deformation_d3(l3_value_global)
end

--************BOUNDARY CONDITIONS*********--
ExtensionEquation_Dirich=DirichletBoundary()
--************INLET BOUNDARY**************--
ExtensionEquation_Dirich:add(0,"k","inlet")
--************WALL BOUNDARY**************--
ExtensionEquation_Dirich:add(0,"k","wall")
--************OUTLET BOUNDARY**************--
ExtensionEquation_Dirich:add(0,"k","outlet")

-- Domain Discretization 
ExtensionEquation_DomainDisc = DomainDiscretization(ExtensionEquation_ApproxSpace)
ExtensionEquation_DomainDisc:add(ExtensionEquation_ElemDisc)
--ExtensionEquation_DomainDisc:add(ExtensionEquation_Dirich)
print("EXTENSION EQUATION: create GridFunctions and GlobalGridFunctionGradientDatas")
k = GridFunction(ExtensionEquation_ApproxSpace);k:set(0.0);k_copy=GridFunction(ExtensionEquation_ApproxSpace);k_copy:set(0.0);
x = GridFunction(ExtensionEquation_ApproxSpace);x:set(nu_ext)
x_limited = GridFunction(ExtensionEquation_ApproxSpace);x_limited:set(nu_ext);--this initial value doesn't matter
r_k = GridFunction(ExtensionEquation_ApproxSpace);r_k:set(0.0)
x_global=GlobalGridFunctionNumberData(x,"k")
A_k = AssembledLinearOperator(ExtensionEquation_DomainDisc)
print("EXTENSION EQUATION: all set:")

print("CREATION OF SOLVERS")
--Creation of Solvers
--deformationState_Solver = util.oo.non_linear_solver(deformationState_DomainDisc,deformationState_ApproxSpace)
deformationState_Solver = util.oo.non_linear_solver_gs(deformationState_DomainDisc,deformationState_ApproxSpace)
--NavierStokes_Solver = util.oo.ns_solver(NavierStokes_DomainDisc,NavierStokes_ApproxSpace)
NavierStokes_Solver = CreateNSSolver(NavierStokes_ApproxSpace,"p","bicgstab",NavierStokes_DomainDisc)

--AdjointFlow_Solver = util.oo.adjoint_ns_solver(AdjointFlow_DomainDisc,AdjointFlow_ApproxSpace)
AdjointFlow_Solver = CreateAdjointFlowSolver(AdjointFlow_ApproxSpace, "h", "bicgstab", AdjointFlow_DomainDisc)

--AdjointDeformation_Solver = util.oo.linear_solver(AdjointDeformation_DomainDisc,AdjointDeformation_ApproxSpace)
AdjointDeformation_Solver = util.oo.linear_solver_gs(AdjointDeformation_DomainDisc,AdjointDeformation_ApproxSpace)

DesignEquation_Solver = util.oo.linear_solver(SurfaceDesignEquation_DomainDisc,SurfaceDesignEquation_ApproxSpace)

ExtensionEquation_Solver = util.oo.linear_solver(ExtensionEquation_DomainDisc,ExtensionEquation_ApproxSpace)

print("DEFORMATION EQUATION SOLVER IS:")
print(deformationState_Solver:config_string())
print("NAVIER STOKES SOLVER IS:")
print(NavierStokes_Solver:config_string())
print("ADJOINT FLOWS SOLVER IS:")
print(AdjointFlow_Solver:config_string())

function ComputeNonLinearSolution(u, domainDisc, solver)
	util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
end

--TODO:PROFILER...no profiling if pdf printing, so set bOutputPlot = false, dos not work in cygwin
function PrintStats() 
	if GetProfilerAvailable() == true then
		
		WriteProfileData("object-optim.pdxml")
		--local scriptPath = "@/afs/math.uni-hamburg.de/users/oa/bay8214/ug4/apps/obstacle_optim/3d_obstacle_optim.lua:"--compu5 Prof. Siebeneborn
		local scriptPath = "@/afs/math.uni-hamburg.de/users/oa/bay8214/ug4/apps/obstacle_optim/3d_obstacle_optim.lua:"--compu5 Jose
		--local scriptPath = "@/home/alfon/ug4/apps/obstacle_optim/3d_obstacle_optim.lua:"--cygwin
		--local scriptPath = "/zhome/academic/HLRS/xbo/xbopinzo/ug4/apps/obstacle_optim/3d_obstacle_optim.lua:"--hawk
		local deformLine = 979 
		local navierLine = 1007
		local adjointFlowLine=1059
		local adjointDefLine=1066
		local controlLine=1071
		local extensionLine = 1079 --search for: extensionLine to get line number
		stats = {
			
			{ "Procs", NumProcs() },
			{ "NumRefs", numRefs },
			    

			{ "DefEq_LS_Init_time_ms", 	GetProfileNode("NewtonPrepareLinSolver", GetProfileNode(scriptPath..(deformLine+2).." ")):get_avg_total_time_ms() },
			{ "DefEq_LS_Apply_time_ms", GetProfileNode("NewtonApplyLinSolver", GetProfileNode(scriptPath..(deformLine+2).." ")):get_avg_total_time_ms() },
			{ "DefEq_LS_Apply_hits", 	GetProfileNode("NewtonApplyLinSolver", GetProfileNode(scriptPath..(deformLine+2).." ")):get_avg_entry_count() },
			{ "DefEq_Newton_Apply_time_ms", GetProfileNode("NewtonSolver_apply", GetProfileNode(scriptPath..(deformLine+2).." ")):get_avg_total_time_ms() },
			{ "DefEq_Newton_Apply_hits", 	GetProfileNode("NewtonSolver_apply", GetProfileNode(scriptPath..(deformLine+2).." ")):get_avg_entry_count() },
			
			{"vSolverDefState_NonLinIterSteps", vSolverDefState_NonLinIterSteps["all"]},
			{"vSolverDefState_LinIterCalls", vSolverDefState_LinIterCalls["all"]},
			{"vSolverDefState_LinIterCalls", vSolverDefState_LinIterCalls["all"]},
			{"vSolverDefState_LinIterAvgSteps", vSolverDefState_LinIterAvgSteps["all"]},
			
			-
			{ "NavEq_LS_Init_time_ms", 	GetProfileNode("NewtonPrepareLinSolver", GetProfileNode(scriptPath..(navierLine+2).." ")):get_avg_total_time_ms() },
			{ "NavEq_LS_Apply_time_ms", GetProfileNode("NewtonApplyLinSolver", GetProfileNode(scriptPath..(navierLine+2).." ")):get_avg_total_time_ms() },
			{ "NavEq_LS_Apply_hits", 	GetProfileNode("NewtonApplyLinSolver", GetProfileNode(scriptPath..(navierLine+2).." ")):get_avg_entry_count() },
			{ "NavEq_Newton_Apply_time_ms", GetProfileNode("NewtonSolver_apply", GetProfileNode(scriptPath..(navierLine+2).." ")):get_avg_total_time_ms() },
			{ "NavEq_Newton_Apply_hits", 	GetProfileNode("NewtonSolver_apply", GetProfileNode(scriptPath..(navierLine+2).." ")):get_avg_entry_count() },
			
			{"vSolverNS_NonLinIterSteps", vSolverNS_NonLinIterSteps["all"]},
			{"vSolverNS_LinIterSteps", vSolverNS_LinIterSteps["all"]},
			{"vSolverNS_LinIterCalls", vSolverNS_LinIterCalls["all"]},
			{"vSolverNS_LinIterAvgSteps", vSolverNS_LinIterAvgSteps["all"]},
			
			
			{ "AdjointFlow_LS_Assemble_time_ms", 	GetProfileNode("assemble_linear", GetProfileNode(scriptPath..(adjointFlowLine+0).." ")):get_avg_total_time_ms() },
			{ "AdjointFlow_LS_Init_time_ms", 		GetProfileNode("LS_InitPrecond", GetProfileNode(scriptPath..(adjointFlowLine+1).." ")):get_avg_total_time_ms() },
			{ "AdjointFlow_LS_Apply_time_ms", 	GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(adjointFlowLine+2).." ")):get_avg_total_time_ms() },
			{ "AdjointFlow_LS_hits", 				GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(adjointFlowLine+2).." ")):get_avg_entry_count() },
			
			{ "vSolverAdjointFlows_LinIterSteps",	vSolverAdjointFlows_LinIterSteps["all"] },
			{ "vSolverAdjointFlows_LinIterAvgRate",	(vSolverAdjointFlows_LinIterAvgRate["all"] / (#vSolverAdjointFlows_LinIterAvgRate+1)) },
			
			
			{ "AdjointDef_LS_Assemble_time_ms", 	GetProfileNode("assemble_linear", GetProfileNode(scriptPath..(adjointDefLine+0).." ")):get_avg_total_time_ms() },
			{ "AdjointDef_LS_Init_time_ms", 		GetProfileNode("LS_InitPrecond", GetProfileNode(scriptPath..(adjointDefLine+1).." ")):get_avg_total_time_ms() },
			{ "AdjointDef_LS_Apply_time_ms", 	GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(adjointDefLine+2).." ")):get_avg_total_time_ms() },
			{ "AdjointDef_LS_hits", 				GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(adjointDefLine+2).." ")):get_avg_entry_count() },
			
			{ "vSolverAdjointDef_LinIterSteps",	vSolverAdjointDef_LinIterSteps["all"] },
			{ "vSolverAdjointDef_LinIterAvgRate",	(vSolverAdjointDef_LinIterAvgRate["all"] / (#vSolverAdjointDef_LinIterAvgRate+1)) },	
			
			{ "Control_LS_Assemble_time_ms", 	GetProfileNode("assemble_linear", GetProfileNode(scriptPath..(controlLine+0).." ")):get_avg_total_time_ms() },
			{ "Control_LS_Init_time_ms", 		GetProfileNode("LS_InitPrecond", GetProfileNode(scriptPath..(controlLine+1).." ")):get_avg_total_time_ms() },
			{ "Control_LS_Apply_time_ms", 	GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(controlLine+2).." ")):get_avg_total_time_ms() },
			{ "Control_LS_hits", 				GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(controlLine+2).." ")):get_avg_entry_count() },

			
			{ "Extension_LS_Assemble_time_ms", 	GetProfileNode("assemble_linear", GetProfileNode(scriptPath..(extensionLine+0).." ")):get_avg_total_time_ms() },
			{ "Extension_LS_Init_time_ms", 		GetProfileNode("LS_InitPrecond", GetProfileNode(scriptPath..(extensionLine+1).." ")):get_avg_total_time_ms() },
			{ "Extension_LS_Apply_time_ms", 	GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(extensionLine+2).." ")):get_avg_total_time_ms() },
			{ "Extension_LS_hits", 				GetProfileNode("LS_ApplyReturnDefect", GetProfileNode(scriptPath..(extensionLine+2).." ")):get_avg_entry_count() },

			
			{ "main", GetProfileNode("root"):get_avg_total_time_ms() },		
		
		}--end stats
		
		if ProcRank() == 0 then
			outdir = util.GetParam("-outdir", "", "output directory")--where do we want to have the measurement file?
			util.printStats(stats)
			util.writeFileStats(stats, outdir.."stats_obstacle_optim.txt")
        else 
            -- sleep 5 seconds
			local clock = os.clock
            local t0 = clock()
            while clock() - t0 <= 5 do end
        end--end procrank==0
		
	end --end if profiler available
end--end print stats

step = 1

vOptimStep = {}

vObjective = {}
vObjective_Drag = {}
vObjective_Control = {}
vObjective_Threshold = {}
vObjective_Volume = {}
vObjective_Barycenter = {}
vObjective_LambdaVolume = {}
vObjective_LambdaBarycenter = {}
vObjective_Extension = {}

vVolume = {}
vBarycenterX = {}
vBarycenterY = {}
vBarycenterZ = {}
vnu_geom = {}
vLambdaBarycenterX = {}
vLambdaBarycenterY = {}
vLambdaBarycenterZ = {}
vLambdaVolume = {}

vControlGradientNorm = {}
vControlNorm = {}

vSolverDefState_NonLinIterSteps = {}
vSolverDefState_NonLinIterSteps["all"] = 0
vSolverDefState_LinIterSteps = {}
vSolverDefState_LinIterSteps["all"] = 0
vSolverDefState_LinIterCalls = {}
vSolverDefState_LinIterCalls["all"] = 0
vSolverDefState_LinIterAvgSteps = {}
vSolverDefState_LinIterAvgSteps["all"] = 0

vSolverNS_NonLinIterSteps = {}
vSolverNS_NonLinIterSteps["all"] = 0
vSolverNS_LinIterSteps = {}
vSolverNS_LinIterSteps["all"] = 0
vSolverNS_LinIterCalls = {}
vSolverNS_LinIterCalls["all"] = 0
vSolverNS_LinIterAvgSteps = {}
vSolverNS_LinIterAvgSteps["all"] = 0

vSolverAdjointFlows_LinIterSteps = {}
vSolverAdjointFlows_LinIterSteps["all"] = 0
vSolverAdjointFlows_LinIterAvgRate = {}
vSolverAdjointFlows_LinIterAvgRate["all"] = 0

vSolverAdjointDef_LinIterSteps = {}
vSolverAdjointDef_LinIterSteps["all"] = 0
vSolverAdjointDef_LinIterAvgRate = {}
vSolverAdjointDef_LinIterAvgRate["all"] = 0

vSolverDesignEq_LinIterSteps = {}
vSolverDesignEq_LinIterSteps["all"] = 0
vSolverDesignEq_LinIterAvgRate = {}
vSolverDesignEq_LinIterAvgRate["all"] = 0

vSolverExtensionEq_LinIterSteps = {}
vSolverExtensionEq_LinIterSteps["all"] = 0
vSolverExtensionEq_LinIterAvgRate = {}
vSolverExtensionEq_LinIterAvgRate["all"] = 0

if bLoadPreviousSteps == true then 
----[[
	print("LOAD VECTOR: FROM CONNECTION VIEWER FILE")
	--LoadVector(d,"DesignEquationSolution.vec")
	--LoadVector(x,"ExtensionEquationSolution.vec")
	ReadFromFile(d,"DesignEquationSolution.vec")
	ReadFromFile(x,"ExtensionEquationSolution.vec")
	print("LOAD VECTOR:SUCCESS")
	--SaveVectorForConnectionViewer(d,"TestD.vec")
	--SaveVectorForConnectionViewer(x,"TestX.vec")
----]]	
end

--G_LBFGS=lbfgs.new()
--G_LBFGS:init(limit_memory,"g","outer","obstacle_surface",SurfaceIntegral)
--G_LBFGS:init_grid_functions(SurfaceDesignEquation_ApproxSpace);


LBFGS=vbfgs.new()
LBFGS:init(limit_memory,"g","outer","obstacle_surface",SurfaceIntegral,"k","outer",nil,VolumeIntegral)
LBFGS:init_grid_functions(SurfaceDesignEquation_ApproxSpace,ExtensionEquation_ApproxSpace);

vtkWriter = VTKOutput()

vtkWriter:print("Grid", w, 0, 0, false)
print_subsets("Surface",w,"obstacle_surface",0,0, false)
print("\nSOLVE PHASE: STARTING ...")
for step=0, numSteps do
print(" +++++++++ Optimization Step " .. step .." +++++++++ ")

vOptimStep[step] = step
boundaryControl_elemDisc:set_boundary_control(d_global)
SurfaceRHS_ElemDisc:set_boundary_control(d_global)
deformationState_elemDisc:set_extension_factor(x_global)
AdjointDeformation_ElemDisc:set_extension_factor(x_global)
ExtensionEquation_ElemDisc:set_extension_factor(x_global)

--print("\NonLinear Solver: solving...")
print("SOLVE PHASE: NON-LINEAR SOLUTION OF THE DEFORMATION STATE EQUATION")
--ComputeNonLinearSolution(w, deformationState_DomainDisc, deformationState_Solver)
deformationState_Solver:init(deformation_Op)--TODO:deformLine
if deformationState_Solver:prepare(w) == false then print ("DeformationStateEquation: Newton solver prepare failed at step "..step.."."); PrintStats(); exit(); end 			
if deformationState_Solver:apply(w) == false then print ("DeformationStateEquation: Newton solver failed at step "..step.."."); PrintStats(); exit(); end
 
--Update the barycenter and volume defects with the newest deformation field iterate w_k
--Volume and barycenter update
vBarycenterDefect = BarycenterDefect(w,deformationNames,"outer",4)
volumeDefect = VolumeDefect(w,ReferenceVolume, "outer", deformationNames,4,false,1,false)
if dim == 3 then 
AdjointDeformation_ElemDisc:set_barycenter_defect(vBarycenterDefect[1],vBarycenterDefect[2],vBarycenterDefect[3])
elseif dim == 2 then
AdjointDeformation_ElemDisc:set_barycenter_defect(vBarycenterDefect[1],vBarycenterDefect[2],0)
end 
AdjointDeformation_ElemDisc:set_volume_defect(volumeDefect)

vVolume[step] = volumeDefect
vBarycenterX[step] = vBarycenterDefect[1]
vBarycenterY[step] = vBarycenterDefect[2]
if dim == 3 then 
vBarycenterZ[step] = vBarycenterDefect[3]
end
gnuplot.write_data("__Geometric_Data.txt", {vOptimStep,vVolume,vBarycenterX,vBarycenterY,vBarycenterZ},false)


print("THE THRESHOLD TERM AFTER DEFORMATION "..step.." IS: "..0.5*beta*Threshold(w,threshold,"outer",deformationNames,4))

print("SOLVE PHASE: NON-LINEAR SOLUTION OF THE NAVIER STOKES PROBLEM")
--ComputeNonLinearSolution(v, NavierStokes_DomainDisc, NavierStokes_Solver)
NavierStokes_Solver:init(navier_Op)----TODO:navierLine
if NavierStokes_Solver:prepare(v) == false then print ("TransformedNavierStokes: Newton solver prepare failed at step "..step.."."); PrintStats(); exit(); end 			
if NavierStokes_Solver:apply(v) == false then print ("TransformedNavierStokes: Newton solver failed at step "..step.."."); PrintStats(); exit(); end 

if step%printInterval==0 then 

if bNavierStokesOutput then
	print("Flows are saved to '" .. flowOutputFile .. "'...")
	--local vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(flowNames, nodalFlows)
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print(flowOutputFile, v, step,step, false)
	vtkWriter:print("Flow_Grid", v,step,step, false) 
	--PRESSURE
	vtkWriter:clear_selection()
	vtkWriter:select_nodal("p", "pressure")
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print("pressureOutputFile", v,step,step, false)
end

if bDeformationSurfaceOutput then
	print("Deformation surface is saved ")
	--local vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(deformationNames, deformationField)
	--vtkWriter:select_all(false)	-- write all variables
	if bPrintSingleStep then
		vtkWriter:print_subsets("DeformationSurface",w,"obstacle_surface",0,0, false)
	else
		vtkWriter:print_subsets("DeformationSurface",w,"obstacle_surface",step,step, false)
	end
end

end



if bDeformationVolumeOutput then
	print("Deformation is saved to '" .. deformationOutputFile .. "'...")
	--local vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(deformationNames, deformationField)
	--vtkWriter:select_all(false)	-- write all variables
	if bPrintSingleStep then
		vtkWriter:print(deformationOutputFile, w, 0, 0, false)
	else
		vtkWriter:print(deformationOutputFile, w, step, step, false)
	end
end

print("SOLVE PHASE: LINEAR SOLUTION OF THE ADJOINT FLOWS PROBLEM")
AdjointFlow_DomainDisc:assemble_linear(A, r_q)--TODO:adjointFlowLine
AdjointFlow_Solver:init(A, q)
AdjointFlow_Solver:apply(q, r_q)



print("SOLVE PHASE: LINEAR SOLUTION OF THE ADJOINT DEFORMATION PROBLEM")
AdjointDeformation_DomainDisc:assemble_linear(A_l, r_l)--TODO: adjointDefLine
AdjointDeformation_Solver:init(A_l, l)
if AdjointDeformation_Solver:apply_return_defect(l, r_l) == false then print ("Adjoint Deformation Equation linear solver failed at step"..step.."."); PrintStats(); exit() ;end

print("SOLVE PHASE: LINEAR SOLUTION OF THE DESIGN EQUATION PROBLEM")
SurfaceDesignEquation_DomainDisc:assemble_linear(A_g, r_g)--TODO:controlLine profiler
DesignEquation_Solver:init(A_g, g)
if DesignEquation_Solver:apply_return_defect(g, r_g) == false then print ("Design Equation linear solver failed at step"..step.."."); PrintStats(); exit() ;end

VecScaleAssign(g,1.0,g)
--VecScaleAssign(g,stepSizeBnd,g)

print("SOLVE PHASE: LINEAR SOLUTION OF THE EXTENSION EQUATION PROBLEM")
ExtensionEquation_DomainDisc:assemble_linear(A_k, r_k)--TODO:extensionLine profiler
ExtensionEquation_Solver:init(A_k, k)
if ExtensionEquation_Solver:apply_return_defect(k, r_k) == false then print ("Extension Equation linear solver failed at step"..step.."."); PrintStats(); exit() ;end

VecScaleAssign(k,1.0,k)
--VecScaleAssign(k,stepSizeExt,k)

print("LBFGS IN OPTIM STEP "..step)

LBFGS:update(d,g,x,k)
VecScaleAdd2(d, 1.0, d, stepSizeBnd, g)
VecScaleAdd2(x, 1.0, x, stepSizeExt, k)	

--G_LBFGS:update_internal(d,g,x,k)
--VecScaleAdd2(d, 1.0, d, 1.0, g)
--VecScaleAdd2(x, 1.0, x, stepSizeExt, k)	

LimitGF(x_limited,x,ext_overwrite,ext_lower);
--Clone into x, because this is connection to the GlobalGridFunctionNumberData x_global
x:assign(x_limited)
if step%printInterval == 0 then

if bAdjointFlowsOutput then
	print("Adjoints is saved to '" .. adjointFlowOutputFile .. "'...")
	--local vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(adjointFlowNames, adjointNodalFlows)
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print(adjointFlowOutputFile, q,step,step, false)
	--PRESSURE
	vtkWriter:clear_selection()
	vtkWriter:select_nodal("h", "adjoint_pressure")
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print("adjointPressureOutputFile", q,step,step, false)
end

if bAdjointDeformationsOutput then
	print("Adjoints is saved to '" .. adjointDeformationOutputFile .. "'...")
	--local vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(adjointDeformationNames, adjointDeformationField)
	--vtkWriter:select_all(false)	-- write all variables
	vtkWriter:print(adjointDeformationOutputFile, l,step,step, false)
end

if bDesignEqOutput then
	print("BoundaryControl Gradient is saved to '" .. designEquationOutputFile .. "'...")
	--vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(designEquationNames, designEquation)
	--vtkWriter:select_all(false)	-- write all variables
	if bPrintSingleStep then
		vtkWriter:print(designEquationOutputFile, d,0,0, false)
	else
		vtkWriter:print(designEquationOutputFile, d,step,step, false)
	end
	
end

if bExtensionEqOutput then
	print("K is saved to '" .. ExtensionEquationOutputFile .. "'...")
	--vtkWriter = VTKOutput()
	vtkWriter:clear_selection()
	vtkWriter:select_nodal(ExtensionEquationNames, ExtensionFactor)
	--vtkWriter:select_all(false)	-- write all variables
	if bPrintSingleStep then
		vtkWriter:print(ExtensionEquationOutputFile, x,0,0, false)
	else
		vtkWriter:print(ExtensionEquationOutputFile, x,step,step, false)
		vtkWriter:print("Extension_Gradient", k,step,step, false)
	end
		
end	

end
	vSolverDefState_NonLinIterSteps[step] = deformationState_Solver:last_num_newton_steps()+1
	vSolverDefState_LinIterSteps[step] = deformationState_Solver:total_linsolver_steps()
	vSolverDefState_LinIterCalls[step] = deformationState_Solver:total_linsolver_calls()
	
	vSolverDefState_NonLinIterSteps["all"] = vSolverDefState_NonLinIterSteps["all"]+deformationState_Solver:last_num_newton_steps()+1
	vSolverDefState_LinIterSteps["all"] = vSolverDefState_LinIterSteps["all"]+deformationState_Solver:total_linsolver_steps()
	vSolverDefState_LinIterCalls["all"] = vSolverDefState_LinIterCalls["all"]+deformationState_Solver:total_linsolver_calls()
	
	vSolverNS_NonLinIterSteps[step] = NavierStokes_Solver:last_num_newton_steps()+1
	vSolverNS_LinIterSteps[step] = NavierStokes_Solver:total_linsolver_steps()
	vSolverNS_LinIterCalls[step] = NavierStokes_Solver:total_linsolver_calls()
	
	vSolverNS_NonLinIterSteps["all"] = vSolverNS_NonLinIterSteps["all"]+NavierStokes_Solver:last_num_newton_steps()+1
	vSolverNS_LinIterSteps["all"] = vSolverNS_LinIterSteps["all"]+NavierStokes_Solver:total_linsolver_steps()
	vSolverNS_LinIterCalls["all"] = vSolverNS_LinIterCalls["all"]+NavierStokes_Solver:total_linsolver_calls()
	
	vSolverAdjointFlows_LinIterSteps[step] =  AdjointFlow_Solver:step()
	vSolverAdjointFlows_LinIterSteps["all"] = vSolverAdjointFlows_LinIterSteps["all"]+ AdjointFlow_Solver:step()
	vSolverAdjointFlows_LinIterAvgRate[step] = AdjointFlow_Solver:convergence_check():avg_rate()

	vSolverAdjointDef_LinIterSteps[step] = AdjointDeformation_Solver:step()
	vSolverAdjointDef_LinIterSteps["all"] = vSolverAdjointDef_LinIterSteps["all"]+ AdjointDeformation_Solver:step()
	vSolverAdjointDef_LinIterAvgRate[step] = AdjointDeformation_Solver:convergence_check():avg_rate()
	
	vSolverDesignEq_LinIterSteps[step] = DesignEquation_Solver:step()
	vSolverDesignEq_LinIterSteps["all"] = vSolverDesignEq_LinIterSteps["all"]+DesignEquation_Solver:step()
	vSolverDesignEq_LinIterAvgRate[step] = DesignEquation_Solver:convergence_check():avg_rate()
	
	vSolverExtensionEq_LinIterSteps[step] = ExtensionEquation_Solver:step()
	vSolverExtensionEq_LinIterSteps["all"] = vSolverExtensionEq_LinIterSteps["all"]+ExtensionEquation_Solver:step()
	vSolverExtensionEq_LinIterAvgRate[step] = ExtensionEquation_Solver:convergence_check():avg_rate()
	
	gnuplot.write_data("__Solver_Stats.txt", {vOptimStep, vSolverDefState_NonLinIterSteps,vSolverDefState_LinIterSteps, vSolverDefState_LinIterCalls,
														  vSolverNS_NonLinIterSteps, vSolverNS_LinIterSteps, vSolverNS_LinIterCalls, 
														  vSolverAdjointFlows_LinIterSteps, vSolverAdjointFlows_LinIterAvgRate, 
														  vSolverAdjointDef_LinIterSteps, vSolverAdjointDef_LinIterAvgRate, 
														  vSolverDesignEq_LinIterSteps, vSolverDesignEq_LinIterAvgRate,
														  vSolverExtensionEq_LinIterSteps, vSolverExtensionEq_LinIterAvgRate
														  },
											 false)
	
	vLambdaBarycenterX[step] = LambdaBarycenter[1]
	vLambdaBarycenterY[step] = LambdaBarycenter[2]
	if dim == 3 then 
	vLambdaBarycenterZ[step] = LambdaBarycenter[3]
	end
	
	vLambdaVolume[step] = LambdaVolume
	vnu_geom[step] = nu_geom
		
	gnuplot.write_data("__Lagrange_Data.txt", {vOptimStep,vLambdaVolume,vLambdaBarycenterX,vLambdaBarycenterY,vLambdaBarycenterZ,vnu_geom},false)
	
	if bOutputPlot then
		gnuplot.plot("_Barycenter.pdf", 
						{
						{ label = "barycenter",  file = "./__Geometric_Data.txt", style = "linespoints", 3, 4 },  
						},
						{	title =	"Barycenter Position", label = {x = "X-coord", y = "Y-coord"}}
					)
		gnuplot.plot("_Volume.pdf", 
					{
						{ label = "volume",  file = "./__Geometric_Data.txt", style = "linespoints", 1, 4 },  
						},
					{	title =	"Volume defect", label = {x = "Step", y = "Defect"}}
					)
	end
	
local drag = 0.5 * visc * Drag(w,v,flowNames,deformationNames,"outer",3)
--TODO:check these terms
--local cont = 0.5 * alpha * Control(d,"g","outer",4)
local cont = 0.5 * alpha * SurfaceIntegral(d,d,"g","obstacle_surface","outer",4)
local thresh = 0.5 * beta * Threshold(w,threshold,"outer",deformationNames,4);
local ext = 0.5 * chi * Extension(x,ext_upper,ext_lower,"outer","k",4);												  
vObjective_Drag[step] = drag
vObjective_Control[step] = cont
vObjective_Threshold[step] = thresh
vObjective_Extension[step] = ext
vObjective_Volume[step]  = 0.5*nu_geom*volumeDefect*volumeDefect
if dim == 2 then 
vObjective_Barycenter[step] = 0.5*nu_geom*
							 (vBarycenterDefect[1]*vBarycenterDefect[1]+vBarycenterDefect[2]*vBarycenterDefect[2])
end 
if dim == 3 then 
vObjective_Barycenter[step] = 0.5*nu_geom*
							 (vBarycenterDefect[1]*vBarycenterDefect[1]+vBarycenterDefect[2]*vBarycenterDefect[2]+vBarycenterDefect[3]*vBarycenterDefect[3])
end 

vObjective_LambdaVolume[step] = LambdaVolume*volumeDefect

if dim == 2 then 
vObjective_LambdaBarycenter[step] = LambdaBarycenter[1]*vBarycenterDefect[1]+LambdaBarycenter[2]*vBarycenterDefect[2]
end
if dim == 3 then 
vObjective_LambdaBarycenter[step] = LambdaBarycenter[1]*vBarycenterDefect[1]+LambdaBarycenter[2]*vBarycenterDefect[2]+LambdaBarycenter[3]*vBarycenterDefect[3]
end
vObjective[step] = drag + cont + thresh + vObjective_Volume[step] + vObjective_Barycenter[step] + ext

					
gnuplot.write_data("__Objective_Function.txt", {vOptimStep,vObjective,vObjective_Drag,vObjective_Control,
											vObjective_Threshold,vObjective_Volume,vObjective_Extension,
											vObjective_Barycenter,vObjective_LambdaVolume,
											vObjective_LambdaBarycenter}, false)
vControlNorm[step]=SurfaceIntegral(d,d,"g","obstacle_surface","outer",4)
vControlGradientNorm[step]=SurfaceIntegral(g,g,"g","obstacle_surface","outer",4)
gnuplot.write_data("__L2_SurfaceNorms.txt",{vOptimStep, vControlNorm, vControlGradientNorm},false)	
	if bOutputPlot then 
		gnuplot.plot("_Objective.pdf", 
						{
							{ label = "all",  file = "./__Objective_Drag.txt", style = "linespoints", 1, 2 }, 
							{ label = "drag",  file = "./__Objective_Drag.txt", style = "linespoints", 1, 3 }, 
							--{ label = "threshold",  file = "./__Objective_Drag.txt", style = "linespoints", 1, 4 },  
						},
						{	title =	"Objective functional", label = {x = "Optimization step", y = "Value objective"}}
				)

	end
	
	if step > 0 and step%updateParam == 0 then
		--Check if subproblem is converged
		if math.abs(vControlNorm[step]-vControlNorm[step-1])/(vControlNorm[step]) < control_tol then
			local tol=math.sqrt(volumeDefect*volumeDefect + vBarycenterDefect[1]*vBarycenterDefect[1] + vBarycenterDefect[2]*vBarycenterDefect[2]+ vBarycenterDefect[3]*vBarycenterDefect[3])
			if tol < geom_constraints_tol then
				print("NO UPDATE OF THE LAGRANGE MULTIPLIERS")
				print("Geometric defect norms is "..tol)
			else
				--Else adjust the Lagrange multipliers
				LambdaVolume = LambdaVolume + multiplier_scaling * nu_geom * volumeDefect
				LambdaBarycenter[1] = LambdaBarycenter[1] + multiplier_scaling * nu_geom * vBarycenterDefect[1]
				LambdaBarycenter[2] = LambdaBarycenter[2] + multiplier_scaling * nu_geom * vBarycenterDefect[2]
		
				LambdaBarycenter[3] = LambdaBarycenter[3] + multiplier_scaling * nu_geom * vBarycenterDefect[3]

				AdjointDeformation_ElemDisc:set_lambda_vol(LambdaVolume)
				
				if dim == 2 then 
				AdjointDeformation_ElemDisc:set_lambda_barycenter(LambdaBarycenter[1],LambdaBarycenter[2],0)	
				end
				if dim == 3 then 
				AdjointDeformation_ElemDisc:set_lambda_barycenter(LambdaBarycenter[1],LambdaBarycenter[2],LambdaBarycenter[3])	
				end		
				--Reset BFGS methods
				--G_LBFGS:reset()
				LBFGS:reset()
			end
			
		end
	end
	

end --step for end

PrintStats(); 


