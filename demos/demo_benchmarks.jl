using BenchmarkTools
using RegLS


vrb = 0
suite = BenchmarkGroup()

solvers = ["ZeroFPR", "FPG", "PG"]
demos   = ["SparseDeconvolution"]

for d in demos
	suite[d] = BenchmarkGroup()
end

include(demos[1]*".jl")
p = SparseDeconvolution.solve_problem!
setup, = SparseDeconvolution.set_up() 
slv_opt = "(verbose = vrb)"

for i in eachindex(solvers)

	solver = eval(parse(solvers[i]*slv_opt))

	suite[demos[1]][solvers[i]] = 
	@benchmarkable p(solver, setup...) setup = (setup = deepcopy($setup); 
						     solver = deepcopy($solver) ) evals = 1
end

results = run(suite, verbose = true)
showall(results)

 
