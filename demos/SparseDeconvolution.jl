
module SparseDeconvolution

using BenchmarkTools
using RegLS
using Convex
using MathProgBase
using Mosek
using AbstractOperators
using RIM
using PyPlot
using PyCall

@pyimport cvxpy as cvx
@pyimport numpy as np

function set_up(;opt="")

	srand(123)
	Fs = 4000        # sampling frequency
	SNR = 15
	env = AcEnv(Fs)

	Lx,Ly,Lz = 4.,5.,3. # room dimensions
	T60 = 0.3           # reverberation time
	geo = CuboidRoom(Lx,Ly,Lz,T60,env) 

	xs = [0.5 0.5 0.5]'            # src pos (in meters)
	xr = [Lx-0.1 Ly-0.3 Lz-0.2]'   # mic pos
	Nh = div(Fs,4)                 # IR samples 
	Nx = div(Fs,2)                 # x  samples 

	t = linspace(0,1/Fs*Nx,Nx)

	h = rim(xs,xr,Nh,geo,env)[:] # Impuse Response

	x = full(sprandn(Nx, 0.05))  # sparse input 

	y0 = conv(x,h)                                      # output signal
	y = y0+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0)) # add noise

	H = Conv(Float64,size(x),h)                              # Abstract Operator
	T = hcat([[zeros(i);h;zeros(Nx-1-i)] for i = 0:Nx-1]...) # Full Matrix
	lambda = 1e-2*vecnorm(H'*y,Inf)

	if opt == "MatrixFree"
		setup = y, H, lambda 
	else
		setup = y, T, lambda 
	end
	return setup, t, x

end


function run_demo()

	setup, t, x = set_up()

	println("Solving Unregularized problem")
	xu = setup[2]\setup[1]

	println("Solving Regularized problem with Full Matrix")
	slv = ZeroFPR()
	x0 = init_variable(x,slv)
	@time x0 = solve_problem!(slv, x0, setup...)
	xm = copy(~x0)

	setup, t, x = set_up(opt = "MatrixFree")

	println("Solving Regularized problem with Abstract Operator")
	~x0 .= 0
	slv = ZeroFPR()
	@time x0 = solve_problem!(slv, x0, setup...)
	x1 = copy(~x0)

	println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
	println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

	return t, x, x1, xu
end

function run_demo_cvx()

	setup, t, x = set_up()

	println("Solving Unregularized problem")
	xu = setup[2]\setup[1]

	#println("Solving Regularized problem with Abstract Operator")
	slv = cvx.SCS
	x0 = init_variable(x,slv)
	@time x0 = solve_problem!(slv, x0, setup...)
	x1 = x0[:value]

	println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
	println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

	return t, x, x1, xu
end

function run_demo_Convex()

	setup, t, x = set_up()

	println("Solving Unregularized problem")
	xu = setup[2]\setup[1]

	println("Solving Regularized problem with Convex")
	slv = MosekSolver()
	x0 = init_variable(x,slv)
	@time x0 = solve_problem!(slv, x0, setup...)
	x1 = x0.value

	println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
	println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

	return t, x, x1, xu
end

init_variable(x,slv::S) where {S <: RegLS.ForwardBackwardSolver} = RegLS.Variable(size(x)...)
init_variable(x,slv::S) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver} = Convex.Variable(size(x)...)
init_variable(x,slv::S) where {S <: AbstractString} = cvx.Variable(size(x)...)

#RegLS Matrix Free
function solve_problem!(slv::S, x0, y, H::A, lambda) where {S <: RegLS.ForwardBackwardSolver, A <: AbstractOperator}
	@minimize ls(H*x0-y)+lambda*norm(x0,1) with slv
	return x0
end

#RegLS non-Matrix Free
function solve_problem!(slv::S, x0, y, T::A, lambda) where {S <: RegLS.ForwardBackwardSolver, A <: AbstractMatrix }
	@minimize ls(T*x0-y)+lambda*norm(x0,1) with slv
	return x0
end

#Convex
function solve_problem!(slv::S, x0, y, T, lambda) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver}
	problem = minimize(0.5*norm(T*x0-y,2)^2+lambda*norm(x0,1)) 
	Convex.solve!(problem,slv)
	return x0
end

#cvxpy
function solve_problem!(slv::S, x0, y, T, lambda) where {S <: AbstractString}
	problem = cvx.Problem(cvx.Minimize(cvx.sum_squares(PyObject(T)*x0-y)*0.5+cvx.norm1(x0)*lambda))
	problem[:solve](solver = slv, verbose = false)
	return x0
end

function benchmark(;verb = 0, samples = 5, seconds = 100)

	suite = BenchmarkGroup()

	solvers = ["ZeroFPR", "FPG", "PG", "cvx.CVXOPT", "cvx.SCS"]
	slv_opt = ["(verbose = $verb)", "(verbose = $verb)", "(verbose = $verb)", "", ""]

	for i in eachindex(solvers)

		setup, t, x = set_up()
		solver = eval(parse(solvers[i]*slv_opt[i]))
		x0 = init_variable(x,solver)

		suite[solvers[i]] = 
		@benchmarkable(solve_problem!(solver, x0, setup...), 
			       setup = (x0 = deepcopy($x0); 
					setup = deepcopy($setup); 
					solver = deepcopy($solver) ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
end

function benchmarkMatrixFree(;verb = 0, samples = 5, seconds = 100)

	suite = BenchmarkGroup()

	solvers = ["ZeroFPR", "FPG", "PG"]
	slv_opt = ["(verbose = $verb)", "(verbose = $verb)", "(verbose = $verb)"]

	for i in eachindex(solvers)

		setup, t, x = set_up(opt = "MatrixFree")
		solver = eval(parse(solvers[i]*slv_opt[i]))
		x0 = init_variable(x,solver)

		suite[solvers[i]] = 
		@benchmarkable(solve_problem!(solver, x0, setup...), 
			       setup = (x0 = deepcopy($x0); 
					setup = deepcopy($setup); 
					solver = deepcopy($solver) ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
end

function show_results(t, x, x1, xu)

	figure()
	plot(t,x )
	plot(t,x1)
	plot(t,xu)

end

end


