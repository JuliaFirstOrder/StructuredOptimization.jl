module LineSpectraEstimation

#TODO change benchmark to JuMP?

using BenchmarkTools
using StructuredOptimization
using Convex
using MathProgBase
using ECOS
using SCS
using AbstractOperators
using PyPlot

function set_up(;opt="")

    srand(39)

    fs = 16e3
    Nt = 2^7 #time samples
    y0 = zeros(Nt)
    SNR = 20

    s = 4 #super-resolution factor
    f_s  = linspace(0,  fs,       s*Nt+1)[1:end-1]     # super resolution frequency axis
    f_s2 = linspace(0,fs/2,div(s*Nt,2)+1)              # super resolution frequency axis (up to Nyquist)
    t  = 0:1/fs:(Nt-1)/fs                              # time axis
    f  = linspace(0,fs,Nt+1)[1:end-1]                  # frequency axis
    f2 = linspace(0,fs/2,div(Nt,2)+1)                  # frequency axis
    K = 14                                             # number of sinusoids
    Nf = div(s*Nt,4)+1
    fk = f_s2[randperm(Nf)[1:K]]                   # sinusoids frequencies
    ak = 0.1*randn(K)+0.7                              # amplitude

    for i in eachindex(fk) y0 .+= ak[i].*sin.(2*π*fk[i].*t) end
    y = y0.+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0))

	xzp = rfft([y;zeros((s-1)*length(y))])
	IDFTm = [exp(im*2*pi*k*n/(s*Nt))  for k =0:s*Nt-1, n=0:s*Nt-1] #Inverse Fourier Matrix
	S = [speye(Nt) spzeros(Nt,(s-1)*Nt)] # selection matrix
	Fm = full(1/(s*Nt).*S*IDFTm)
    alpha = 0.001
	lambda_max_m = norm(Fm'*y, Inf)
	lambda_m = alpha*lambda_max_m

    F = 1/(s*Nt)*DFT(s*Nt)'[1:Nt] # Abstract Operator 
	lambda_max = norm(F'*y, Inf)
	lambda = alpha*lambda_max

	Fc = [real(Fm) -imag(Fm); imag(Fm) real(Fm)] 
	# needed in cvxpy (currently does not support complex variables)

	if opt == "MatrixFree"
		setup = K,  F, Fc, lambda, lambda_m
	else
		setup = K, Fm, Fc, lambda, lambda_m
	end
	return setup, t, f, f_s2, fk, ak, s, Nt, fs, xzp, y

end

function run_demo()

    tol = 1e-8
    verb = 1
    solver = PANOC

	setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up(opt = "MatrixFree")
	slv = solver(tol = tol, verbose = verb)
	x = init_variable(s*Nt,slv)

	println("Solving LASSO (Abstract Operator)")
	@time solve_problem!(slv, x, y, setup...)
	x1 = copy(~x)

	println("Refine solution by solving non-convex problem (Abstract Operator)")
	@time solve_problem_ncvx!(slv, x, y, setup...)
	x0 = copy(~x)

	slv = solver(tol = tol, verbose = verb)
	setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up()
	x0m = init_variable(s*Nt,slv)

	println("Solving LASSO (Matrix Operator)")
	@time solve_problem!(slv, x0m, y, setup...)

	println("Refine solution by solving non-convex problem (Matrix Operator)")
	@time solve_problem_ncvx!(slv, x0m, y, setup...)

	x1 = x1[1:div(s*Nt,2)+1]
	x0 = x0[1:div(s*Nt,2)+1]
	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, x0
end

function run_demo_Convex()
	slv = ECOSSolver()
	setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up()
	x0m = init_variable(s*Nt,slv)

	println("Solving LASSO with Convex.jl (Matrix Operator)")
	@time x0m, = solve_problem!(slv, x0m, y, setup...)
	
	x1 = x0m.value
	x1 = x1[1:div(s*Nt,2)+1]
	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, zeros(x1)
end

init_variable(N,slv::S) where {S <: StructuredOptimization.ForwardBackwardSolver} = StructuredOptimization.Variable(Complex{Float64}, N)
init_variable(N,slv::S) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver} = Convex.ComplexVariable(N)
init_variable(N,slv::S) where {S <: AbstractString} = cvx.Variable(N), cvx.Variable(N)

#StructuredOptimization Matrix Free
function solve_problem!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: StructuredOptimization.ForwardBackwardSolver, A <: AbstractOperator}
    it, = @minimize ls(F*x0-y)+lambda*norm(x0,1) with slv
	return x0, it
end

function solve_problem_ncvx!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: StructuredOptimization.ForwardBackwardSolver, A <: AbstractOperator}
	it, = @minimize ls(F*x0-y) st norm(x0,0) <= 2*K with slv
	return x0, it
end

#StructuredOptimization non-Matrix Free
function solve_problem!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: StructuredOptimization.ForwardBackwardSolver, A <: AbstractMatrix}
    it, = @minimize ls(F*x0-complex(y))+lambda_m*norm(x0,1) with slv
	return x0, it
end

function solve_problem_ncvx!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: StructuredOptimization.ForwardBackwardSolver, A <: AbstractMatrix}
    it, = @minimize ls(F*x0-complex(y)) st norm(x0,0) <= 2*K with slv
	return x0, it
end

#Convex 
function solve_problem!(slv::S, x0, y, K, F, Fc, lambda, lambda_m) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver}
	problem = minimize(0.5*norm(F*x0-y,2)^2+lambda_m*norm(x0,1)) 
	Convex.solve!(problem,slv)
	return x0, 0
end

function benchmark(;verb = 0, samples = 5, seconds = 100, tol = 1e-8, maxit = 20000)

	suite = BenchmarkGroup()

	solvers = ["PANOC", "ZeroFPR", "FPG", "PG",]
	slv_opt = ["(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)",]

	its = Dict([(sol,0.) for sol in solvers])
	for i in eachindex(solvers)

		setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up()
		solver = eval(parse(solvers[i]*slv_opt[i]))
		x0 = init_variable(Nt*s,solver)

		suite[solvers[i]] = 
		@benchmarkable((x0,it) = solve_problem!(solver, x0, y, setup...), 
			       setup = (
					it = 0;
					x0 = deepcopy($x0); 
					setup = deepcopy($setup); 
					y = $y; 
					solver = deepcopy($solver) ), 
			       teardown = (
					  $its[$solvers[$i]] = it;
					  ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
	println("LineSpectraEstimation its")
	println(its)
	return results
end

function benchmarkMatrixFree(;verb = 0, samples = 5, seconds = 100, tol = 1e-8, maxit = 20000)

	suite = BenchmarkGroup()

	solvers = ["PANOC", "ZeroFPR", "FPG", "PG"]
	slv_opt = ["(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)"]

	its = Dict([(sol,0.) for sol in solvers])
	for i in eachindex(solvers)

		setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up(opt = "MatrixFree")
		solver = eval(parse(solvers[i]*slv_opt[i]))
		x0 = init_variable(Nt*s,solver)

		suite[solvers[i]] = 
		@benchmarkable((x0,it) = solve_problem!(solver, x0, y, setup...), 
			       setup = (
					it = 0;
					x0 = deepcopy($x0); 
					setup = deepcopy($setup); 
					y = $y; 
					solver = deepcopy($solver) ), 
			       teardown = (
					  $its[$solvers[$i]] = it;
					  ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
	println("LineSpectraEstimation (MatrixFree) its")
	println(its)
	return results
end

function show_results(t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, x0)

	figure()
	plot(fs,abs.(xzp./Nt ), label = "dft zero pad.")
	plot(f,        abs.(fft(y)./Nt)       , label = "dft")
	plot(fk,       abs.(ak)/2     , "r*", label = "true amp.")
    plot(fs,abs.(x1./(s*Nt)), "k*", label = "LASSO")
    plot(fs,abs.(x0./(s*Nt)), "go", label = "IndBallL0")
	xlim([0;Fs/4])
	legend()

end

end
