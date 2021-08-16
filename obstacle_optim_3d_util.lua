  
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


util.oo = util.oo or {}

--change to baseSolver=LU(),

function util.oo.non_linear_solver(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			="jac",--ComponentGaussSeidel(0.1, {"p"}, {0}, {1}),
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = LU(),--"lu",   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 3000,		-- number of iterations
			absolute	= 1e-14,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-2,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	NonLinSolverDesc = 
	{
	  type    = "newton",
	  convCheck = {
			type		= "standard",
			iterations	= 20,		-- number of iterations
			absolute	= 1e-12,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-8,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   },
	  lineSearch  = 
	  {
	  type      = "standard",
		maxSteps      = 50,
		lambdaStart   = 1,
		lambdaReduce  = 0.9,
		acceptBest    = true,
		checkAll      = false,
		verbose		  = false
	  },
	  linSolver	= LinSolverDesc
	}
	non_LinSolver=util.solver.CreateSolver(NonLinSolverDesc)
	return non_LinSolver
end

function util.oo.non_linear_solver_gs(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			="gs",--ComponentGaussSeidel(0.1, {"p"}, {0}, {1}),
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = LU(),--"lu",   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 3000,		-- number of iterations
			absolute	= 1e-14,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-2,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	NonLinSolverDesc = 
	{
	  type    = "newton",
	  convCheck = {
			type		= "standard",
			iterations	= 20,		-- number of iterations
			absolute	= 1e-12,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-8,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   },
	  lineSearch  = 
	  {
	  type      = "standard",
		maxSteps      = 50,
		lambdaStart   = 1,
		lambdaReduce  = 0.9,
		acceptBest    = true,
		checkAll      = false,
		verbose		  = false
	  },
	  linSolver	= LinSolverDesc
	}
	non_LinSolver=util.solver.CreateSolver(NonLinSolverDesc)
	return non_LinSolver
end

function util.oo.linear_solver(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			="jac",
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = LU(),--"lu",   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 2000,		-- number of iterations
			absolute	= 1e-16,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-16,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	LinSolver=util.solver.CreateSolver(LinSolverDesc)
	return LinSolver
end

function util.oo.linear_solver_gs(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			="gs",
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = LU(),--"lu",   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 2000,		-- number of iterations
			absolute	= 1e-16,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-16,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	LinSolver=util.solver.CreateSolver(LinSolverDesc)
	return LinSolver
end

function util.oo.ns_solver(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			=ComponentGaussSeidel(0.1, {"p"}, {0}, {1}),
			--smoother			=BlockGaussSeidel(1),
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = LU(),--"LU",   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 20000,		-- number of iterations
			absolute	= 1e-10,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-8,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	NonLinSolverDesc = 
	{
	  type    = "newton",
	  convCheck = {
			type		= "standard",
			iterations	= 50,		-- number of iterations
			absolute	= 1e-8,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-6,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   },
	  lineSearch  = 
	  {
	  type      = "standard",
		maxSteps      = 50,
		lambdaStart   = 1,
		lambdaReduce  = 0.9,
		acceptBest    = true,
		checkAll      = false,
		verbose		  = false
	  },
	  linSolver	= LinSolverDesc
	}
	non_LinSolver=util.solver.CreateSolver(NonLinSolverDesc)
	return non_LinSolver
end

function util.oo.adjoint_ns_solver(domainDisc, approxSpace)
	LinSolverDesc = 
	{
		type = "bicgstab",   
		precond=
		{
			type                = "gmg",
			smoother			=ComponentGaussSeidel(0.1, {"h"}, {0}, {1}),
			--smoother			=BlockGaussSeidel(1),
			adaptive            = false,
			approxSpace           = approxSpace,
			baseLevel             = 0,
			gatheredBaseSolverIfAmbiguous = false,
			baseSolver            = LU(),   -- any solver listed in the 'Solvers' section
			--baseSolver            = LU(),
			cycle               = "V",
			discretization          = domainDisc,    -- only necessary if the underlying matrix is not of type AssembledLinearOperator
			preSmooth             = 3,
			postSmooth            = 3,
			rap               = true,
			transfer            = "std",  -- any transfer listed in the 'Transfers' section
			debug               = false,
			mgStats             = nil , -- any mgStats listed in the 'MGStats' section
		},
	   convCheck = {
			type		= "standard",
			iterations	= 20000,		-- number of iterations
			absolute	= 1e-10,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
			reduction	= 1e-10,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
			verbose		= true,		-- print convergence rates if true
	   }
	  --approxSpace = approxSpace
	}
	LinSolver=util.solver.CreateSolver(LinSolverDesc)
	return LinSolver
end

function util.oo.gmg(approxSpace, smoother, numPreSmooth, numPostSmooth,
						 cycle_str, baseSolver, baseLev, bRAP)

	local gmg = GeometricMultiGrid(approxSpace)
	
	gmg:set_base_level(baseLev)
	gmg:set_base_solver(baseSolver)
	gmg:set_gathered_base_solver_if_ambiguous(true)
	gmg:set_smoother(smoother)
	gmg:set_cycle_type(cycle)
	gmg:set_num_presmooth(numPreSmooth)
	gmg:set_num_postsmooth(numPostSmooth)
	gmg:set_rap(bRAP)
	
	return gmg

end
