
module SparseDeconvolution

using BenchmarkTools
using StructuredOptimization
using JuMP, MathProgBase, SCS, ECOS
using AbstractOperators
using RIM
using PyPlot


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

    x_test = randn(Nx)
    @assert norm(T*x_test-H*x_test) < 1e-8

	if opt == "MatrixFree"
		setup = y, H, lambda, Nx 
	elseif opt == "JuMP"
		M = Model()
		@variables M begin
			x0[1:Nx]
			tt[1:Nx]
			w
		end
		@objective(M,Min,[0.5;lambda*ones(Nx)]'*[w;tt])
		@constraint(M, soc, norm( [1-w;2*(T*x0-y)] ) <= 1+w)
		@constraint(M,  x0 .<= tt)
		@constraint(M, -tt .<= x0)
		setup = y, T, M, x0
	else
		setup = y, T, lambda, Nx 
	end
	return setup, t, x

end


function run_demo()

	setup, t, x = set_up()
    tol = 1e-6
    verb = 1
    solver = PANOC

	println("Solving Unregularized problem")
	xu = setup[2]\setup[1]

	println("Solving Regularized problem with Full Matrix")
	slv = solver(tol = tol, verbose = verb)
	@time x0, = solve_problem(slv, setup...)
	xm = copy(~x0)

	setup, t, x = set_up(opt = "MatrixFree")

	println("Solving Regularized problem with Abstract Operator")
	slv = solver(tol = tol, verbose = verb)
	@time x0, = solve_problem(slv, setup...)
	x1 = copy(~x0)

	println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
	println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

	return t, x, x1, xu
end

function run_demo_JuMP()

	setup, t, x = set_up(opt = "JuMP")

	println("Solving Unregularized problem")
	xu = setup[2]\setup[1]

	println("Solving Regularized problem with conic solver")
	slv = SCSSolver()
	@time x1, = solve_problem(slv, setup...)

	println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
	println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

	return t, x, x1, xu
end

#StructuredOptimization Matrix Free
function solve_problem(slv::S, y, H::A, lambda, Nx) where {S <: StructuredOptimization.ForwardBackwardSolver, A <: AbstractOperator}
	x0 = StructuredOptimization.Variable(Nx) 
	it, = @minimize ls(H*x0-y)+lambda*norm(x0,1) with slv
	return x0, it
end

#StructuredOptimization non-Matrix Free
function solve_problem(slv::S, y, T::A, lambda, Nx) where {S <: StructuredOptimization.ForwardBackwardSolver, A <: AbstractMatrix }
	x0 = StructuredOptimization.Variable(Nx) 
	it, = @minimize ls(T*x0-y)+lambda*norm(x0,1) with slv
	return x0, it
end

#JuMP non-Matrix Free
function solve_problem(slv::S, y, T, M, x0) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver}
	JuMP.setsolver(M, slv)
	JuMP.solve(M)
	return getvalue(x0), 0
end

function benchmark(;verb = 0, samples = 5, seconds = 100, tol = 1e-6)

	suite = BenchmarkGroup()

	opt     = [
           #"JuMP",
		   #"JuMP",
		   "",
		   "",
		   "",
		   "",
          ]
	solvers = [
		   #"ECOSSolver",
		   #"SCSSolver",
		   "PANOC",
		   "ZeroFPR", 
		   "FPG", 
		   "PG"]
	slv_opt = [
		   "(verbose = $verb, tol = $tol)", 
		   "(verbose = $verb, tol = $tol)", 
		   "(verbose = $verb, tol = $tol)", 
		   "(verbose = $verb, tol = $tol)", 
		   "(verbose = $verb, tol = $tol)", 
		   "(verbose = $verb, tol = $tol)"]

	its = Dict([(sol,0.) for sol in solvers])
	for i in eachindex(solvers)

		setup, t, x = set_up(opt = opt[i])
		solver = eval(parse(solvers[i]*slv_opt[i]))

		suite[solvers[i]] = 
		@benchmarkable((x0,it) = solve_problem(solver, setup...), 
			       setup = ( 
					it = 0;
					setup = deepcopy($setup); 
					solver = deepcopy($solver) ), 
			       teardown = (
					  $its[$solvers[$i]] = it;
					  ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
	println("SparseDeconvolution its")
	println(its)
	return results
end

function benchmarkMatrixFree(;verb = 0, samples = 5, seconds = 100, tol = 1e-6)

	suite = BenchmarkGroup()

	solvers = ["PANOC", "ZeroFPR", "FPG", "PG"]
	slv_opt = ["(verbose = $verb, tol = $tol)", 
               "(verbose = $verb, tol = $tol)", 
               "(verbose = $verb, tol = $tol)", 
               "(verbose = $verb, tol = $tol)"]

	its = Dict([(sol,0.) for sol in solvers])
	for i in eachindex(solvers)

		setup, t, x = set_up(opt = "MatrixFree")
		solver = eval(parse(solvers[i]*slv_opt[i]))

		suite[solvers[i]] = 
		@benchmarkable((x0,it) = solve_problem(solver, setup...), 
			       setup = ( 
					it = 0;
					setup = deepcopy($setup); 
					solver = deepcopy($solver) ), 
			       teardown = (
					  $its[$solvers[$i]] = it;
					  ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
	println("SparseDeconvolution Matrix Free its")
	println(its)
	return results
end

function show_results(t, x, x1, xu)

	figure()
	plot(t,x )
	plot(t,x1)
	plot(t,xu)

end

end


