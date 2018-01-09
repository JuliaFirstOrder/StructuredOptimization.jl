module LineSpectraEstimation

using BenchmarkTools
using RegLS
using Convex
using MathProgBase
using Mosek
using AbstractOperators
using PyPlot
using PyCall

@pyimport cvxpy as cvx
@pyimport numpy as np

function set_up(;opt="")

	srand(17)

	Fs = 16e3
	Nt = 2^8 #time samples
	y0 = zeros(Nt)
	SNR = 10

	s = 6 #super-resolution factor
	fs = linspace(0,Fs/2,div(s*Nt,2)+1)            # super resolution frequency axis (up to Nyquist)
	t  = 0:1/Fs:(Nt-1)/Fs                          # time axis
	f  = linspace(0,Fs,Nt+1)[1:end-1]              # frequency axis
	K = 14                                         # number of sinusoids
	fk = fs[randperm(div(s*Nt,2)+1)[1:K]]          #sinusoids frequencies
	ak = 0.1*randn(K)+0.7                          # amplitude

	for i in eachindex(fk) y0 .+= ak[i].*sin.(2*π*fk[i].*t) end
	y = y0.+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0))

	xzp = rfft([y;zeros((s-1)*length(y))])
	IDFTm = [exp(im*2*pi*k*n/(s*Nt))  for k =0:s*Nt-1, n=0:s*Nt-1] #Inverse Fourier Matrix
	S = [speye(Nt) spzeros(Nt,(s-1)*Nt)] # selection matrix
	Fm = full(S*IDFTm)
	lambda_max_m = norm(Fm'*y, Inf)
	lambda_m = 0.06*lambda_max_m

	F = (s*Nt*IRDFT((div(s*Nt,2)+1,),s*Nt))[1:Nt] # Abstract Operator 
	lambda_max = norm(F'*y, Inf)
	lambda = 0.06*lambda_max

	Fc = [real(Fm) -imag(Fm); imag(Fm) real(Fm)] 
	# needed in cvxpy (currently does not support complex variables)

	if opt == "MatrixFree"
		setup = K,  F, Fc, lambda, lambda_m
	else
		setup = K, Fm, Fc, lambda, lambda_m
	end
	return setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y

end

function run_demo()
	slv = ZeroFPR()
	setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up(opt = "MatrixFree")
	x = init_variable(div(s*Nt,2)+1,slv)

	println("Solving LASSO (Abstract Operator)")
	@time solve_problem!(slv, x, y, setup...)
	x1 = copy(~x)

	println("Refine solution by solving non-convex problem (Abstract Operator)")
	@time solve_problem_ncvx!(slv, x, y, setup...)
	x0 = copy(~x)

	slv = ZeroFPR()
	setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up()
	x0m = init_variable(s*Nt,slv)

	println("Solving LASSO (Matrix Operator)")
	@time solve_problem!(slv, x0m, y, setup...)

	println("Refine solution by solving non-convex problem (Matrix Operator)")
	@time solve_problem_ncvx!(slv, x0m, y, setup...)

	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, x0
end

function run_demo_cvx()

	slv = cvx.MOSEK
	setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up()
	x0 = init_variable(s*Nt,slv)

	@time x0 = solve_problem!(slv, x0, y, setup...)
	x1 = x0[1][:value]+im*x0[2][:value]
	x1 = x1[1:div(s*Nt,2)+1]

	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, zeros(x1)
end

function run_demo_Convex()
	slv = MosekSolver()
	setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up()
	x0m = init_variable(s*Nt,slv)

	println("Solving LASSO with Convex.jl (Matrix Operator)")
	@time x0m = solve_problem!(slv, x0m, y, setup...)
	
	x1 = x0m.value
	x1 = x1[1:div(s*Nt,2)+1]
	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, zeros(x1)
end

init_variable(N,slv::S) where {S <: RegLS.ForwardBackwardSolver} = RegLS.Variable(Complex{Float64}, N)
init_variable(N,slv::S) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver} = Convex.ComplexVariable(N)
init_variable(N,slv::S) where {S <: AbstractString} = cvx.Variable(N), cvx.Variable(N)

#RegLS Matrix Free
function solve_problem!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: RegLS.ForwardBackwardSolver, A <: AbstractOperator}
	@minimize ls(F*x0-y)+lambda*norm(x0,1) with slv
	return x0
end

function solve_problem_ncvx!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: RegLS.ForwardBackwardSolver, A <: AbstractOperator}
	@minimize ls(F*x0-y) st norm(x0,0) <= K with slv
	return x0
end

#RegLS non-Matrix Free
function solve_problem!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: RegLS.ForwardBackwardSolver, A <: AbstractMatrix}
	@minimize ls(F*x0-y)+lambda_m*norm(x0,1) with slv
	return x0
end

function solve_problem_ncvx!(slv::S, x0, y, K, F::A, Fc, lambda, lambda_m) where {S <: RegLS.ForwardBackwardSolver, A <: AbstractMatrix}
	@minimize ls(F*x0-y) st norm(x0,0) <= 2*K with slv
	return x0
end

#Convex 
function solve_problem!(slv::S, x0, y, K, F, Fc, lambda, lambda_m) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver}
	problem = minimize(0.5*norm(F*x0-y,2)^2+lambda_m*norm(x0,1)) 
	Convex.solve!(problem,slv)
	return x0
end

#cvxpy 
function solve_problem!(slv::S, x0, y, K, F, Fc, lambda, lambda_m) where {S <: AbstractString}
	xr, xi = x0[1],x0[2]
	reg = cvx.norm(cvx.vstack(xr[1],xi[1])) 
	for i = 2:xr[:size][1]
		reg += cvx.norm(cvx.vstack(xr[i],xi[i]))
	end

	problem = cvx.Problem(
		  cvx.Minimize(cvx.sum_squares(PyObject(Fc)*cvx.vstack(xr,xi)-[real(y);imag(y)])*0.5
			       +reg*lambda_m
			       ))
	problem[:solve](solver = slv, verbose = false)
	return x0
end

function benchmark(;verb = 0, samples = 5, seconds = 100)

	suite = BenchmarkGroup()

	solvers = ["ZeroFPR", "FPG", "PG", "cvx.MOSEK", "cvx.SCS"]
	slv_opt = ["(verbose = $verb)", "(verbose = $verb)", "(verbose = $verb)", "", ""]

	for i in eachindex(solvers)

		setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up()
		solver = eval(parse(solvers[i]*slv_opt[i]))
		x0 = init_variable(Nt*s,solver)

		suite[solvers[i]] = 
		@benchmarkable(solve_problem!(solver, x0, y, setup...), 
			       setup = (x0 = deepcopy($x0); 
					setup = deepcopy($setup); 
					y = $y; 
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

		setup, t, f, fs, fk, ak, s, Nt, Fs, xzp, y = set_up(opt = "MatrixFree")
		solver = eval(parse(solvers[i]*slv_opt[i]))
		x0 = init_variable(div(Nt*s,2)+1,solver)

		suite[solvers[i]] = 
		@benchmarkable(solve_problem!(solver, x0, y, setup...), 
			       setup = (x0 = deepcopy($x0); 
					setup = deepcopy($setup); 
					y = $y; 
					solver = deepcopy($solver) ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
end

function show_results(t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, x0)

	figure()
	plot(fs,abs.(xzp./Nt ), label = "dft zero pad.")
	plot(f,        abs.(fft(y)./Nt)       , label = "dft")
	plot(fk,       abs.(ak)/2     , "r*", label = "true amp.")
	plot(fs,abs.(x1), "k*", label = "LASSO")
	plot(fs,abs.(x0), "go", label = "IndBallL0")
	xlim([0;Fs/2])
	legend()

end

end
