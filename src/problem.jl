include("problems/extract.jl")
include("problems/split.jl")

export problem, minimize

function problem(terms...)
	cf = vcat(terms...)
end

minimize(cf::Vararg{Term}; solver::Solver=default_solver()) =
	solve(problem(cf...), solver)
