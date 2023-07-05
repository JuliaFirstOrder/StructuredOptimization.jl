module TotalVariation

using BenchmarkTools
using StructuredOptimization
using AbstractOperators
using Images
using ImageView
using TestImages
using Random

function set_up()
	Random.seed!(1234)
	img = testimage("cameraman")

	X = convert(Array{Float64},img) # convert image to array
	Xt = X.+sqrt(0.006*norm(X,Inf))*randn(size(X)) # add noise
	Xt[Xt .< 0] .= 0. #make sure pixels are in range
	Xt[Xt .> 1] .= 1.

	V = Variation(size(X)) # linear mapping operator
	lambda = 0.07          # level of regularization

	Y = Variable(size(V,1)...) # dual variables
	return V, Y, Xt, X, lambda
end

function run_demo()
	Lf  = 8                # Lipschitz constant
	slv = PANOC(tol = 1e-3, gamma = 1/Lf, adaptive = false)
	setup = set_up()
	@time solve_problem!(slv,setup...)
	return setup
end

function solve_problem!(slv,V, Y, Xt, X, lambda)
	it  = @minimize ls(-V'*Y+Xt)+conj(lambda*norm(Y,2,1,2)) with slv
	return it
end

function benchmark(;verb = 0, samples = 5, seconds = 100, tol = 1e-3, maxit = 50000 )

	suite = BenchmarkGroup()

	solvers = ["ZeroFPR",
               "PANOC",
               "ForwardBackward"]
	slv_opt = ["(verbose = $verb, tol = $tol, gamma = 1/8, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, gamma = 1/8, maxit = $maxit)",
               "(verbose = $verb, tol = $tol, gamma = 1/8, maxit = $maxit)"]

	its = Dict([(sol,0.) for sol in solvers])
	for i in eachindex(solvers)

		setup = set_up()
		solver = eval(Meta.parse(solvers[i]*slv_opt[i]))

		suite[solvers[i]] = 
		@benchmarkable(
			it = solve_problem!(solver, setup...), 
			setup = ( 
				it = 0;
				setup = deepcopy($setup); 
				solver = deepcopy($solver) 
			), 
			teardown = (
				$its[$solvers[$i]] = it[2];
			), 
			evals = 1, samples = samples, seconds = seconds)
	end

	println("Starting run")
	results = run(suite, verbose = (verb != 0))
	println("TotalVariation its")
	println(its)
	return results
end

function show_results(V, Y, Xt, X, lambda)
	Xd = Xt-V'*(~Y)
	imshow(X)
	imshow(Xt)
	imshow(Xd)
end

end
