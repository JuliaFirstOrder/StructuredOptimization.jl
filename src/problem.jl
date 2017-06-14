include("problems/extract.jl")
include("problems/split.jl")
include("problems/solve.jl")

export problem, minimize

function problem(terms::Vararg)
	cf = ()
	for i = 1:length(terms)
		typeof(terms[i]) <: Tuple ? cf = (cf...,terms[i]...) : cf = (cf...,terms[i])
	end
	return cf
end

minimize(cf::Vararg{Term}; solver::Solver=default_solver()) =
	solve(problem(cf...), solver)
