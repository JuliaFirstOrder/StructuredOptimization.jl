export problem, minimize

function problem(terms...)
	cf = vcat(terms...)
	return cf
end

minimize(cf::Vararg{Term}; solver::Solver=default_solver()) =
	solve(problem(cf...), solver)
