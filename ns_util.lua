  
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

-- Create util namespace
util = util or {}


util.ns = util.ns or {}

function util.solver.create(sol, precond)

	local solver = nil
	
	if 		sol == "ls" 		then 
		solver = LinearSolver();
		if precond ~= nil then solver:set_preconditioner(precond) end
	elseif 	sol == "bicgstab" 	then 
		solver = BiCGStab();
		solver:set_min_orthogonality(1e-20)
		--solver:set_restart(30)
		if precond ~= nil then solver:set_preconditioner(precond) end
	elseif 	sol == "gmres" 	then 
		solver = GMRES(10);
		if precond ~= nil then solver:set_preconditioner(precond) end
	elseif 	sol == "cg" 		then 
		solver = CG();
		if precond ~= nil then solver:set_preconditioner(precond) end
	elseif  sol == "schur" then 

		local skeletonSolver = BiCGStab()
		skeletonSolver:set_preconditioner(ILU())
		--skeletonSolver:set_min_orthogonality(1e-15)
		--skeletonSolver:set_restart(30)		
		skeletonSolver:set_convergence_check(ConvCheck(10000, 1e-12, 1e-2, true))	

		skeletonSolver = AgglomeratingSolver(LU())

		local schur = SchurComplement()
		schur:set_dirichlet_solver(LU())
		schur:set_skeleton_solver(SchurInverseWithFullMatrix(skeletonSolver))
	
		solver = LinearSolver()
		solver:set_preconditioner(schur)
			
	elseif sol == "lu" then
		solver = LU()
	else
		print("Solver not found."); exit();
	end	
	
	return solver
end

util.gmg = util.gmg or {}

function util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
						 cycle_string, baseSolver, baseLev, bRAP,discretization)

	local gmg = GeometricMultiGrid(approxSpace)
	
	gmg:set_base_level(baseLev)
	gmg:set_base_solver(baseSolver)
	gmg:set_gathered_base_solver_if_ambiguous(false)
	gmg:set_smoother(smoother)
	gmg:set_cycle_type(cycle_string)
	gmg:set_num_presmooth(numPreSmooth)
	gmg:set_num_postsmooth(numPostSmooth)
	gmg:set_rap(bRAP)
	gmg:set_discretization(discretization)
	return gmg
end