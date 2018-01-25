module TotalVariation

using BenchmarkTools
using RegLS
using AbstractOperators
using Images
using ImageView
using TestImages

function set_up()
	srand(123)
	img = testimage("cameraman")

	X = convert(Array{Float64},img) # convert image to array
	Xt = X.+sqrt(0.006*vecnorm(X,Inf))*randn(size(X)) # add noise
	Xt[Xt .< 0] .= 0. #make sure pixels are in range
	Xt[Xt .> 1] .= 1.

	V = Variation(size(X)) # linear mapping operator
	lambda = 0.07          # level of regularization

	Y = Variable(size(V,1)...) # dual variables
	return V, Y, Xt, X, lambda
end

function run_demo()
	Lf  = 8                # Lipschitz constant
	slv = ZeroFPR(tol = 1e-3, gamma = 1/Lf, adaptive = false)
	setup = set_up()
	@time solve_problem!(slv,setup...)
	return setup
end

function solve_problem!(slv,V, Y, Xt, X, lambda)
	@minimize ls(V'*Y-Xt)+conj(lambda*norm(Y,2,1,2)) with slv
end

function benchmark(;verb = 0, samples = 5, seconds = 100)

	suite = BenchmarkGroup()

	tol = 1e-3
	solvers = ["ZeroFPR",
		   "FPG",
		   "PG"]
	slv_opt = ["(verbose = $verb, tol = $tol, gamma = 1/8, maxit = 50000)", 
		   "(verbose = $verb, tol = $tol, gamma = 1/8, maxit = 50000)",
		   "(verbose = $verb, tol = $tol, gamma = 1/8, maxit = 50000)"]

	for i in eachindex(solvers)

		setup = set_up()
		solver = eval(parse(solvers[i]*slv_opt[i]))

		suite[solvers[i]] = 
		@benchmarkable(solve_problem!(solver, setup...), 
			       setup = ( 
					setup = deepcopy($setup); 
					solver = deepcopy($solver) ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
end

function show_results(V, Y, Xt, X, lambda)
	Xd = -(V'*(~Y)-Xt)
	imshow(X)
	imshow(Xt)
	imshow(Xd)
end

end
