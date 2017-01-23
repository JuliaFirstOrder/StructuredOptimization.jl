__precompile__()

module RegLS

using ProximalOperators

include("solvers.jl")

export solve, solve!

solve(A, b::AbstractArray, g::ProximableFunction) =
	solve(A, b::AbstractArray, g::ProximableFunction, zeros(eltype(b), size(A,2)))

function solve(A, b::AbstractArray, g::ProximableFunction, args...)
	L = x -> A*x
	Ladj = x -> A'*x
	solve(L, Ladj, b, g, args...)
end

solve(L::Function, Ladj::Function, b::AbstractArray, g::ProximableFunction, x::Array) =
  solve(L, Ladj, b, g, x, ZeroFPR())

function solve(L::Function, Ladj::Function, b::AbstractArray, g::ProximableFunction, x0::AbstractArray, args...)
	x = deepcopy(x0) # copy initial conditions
	x, slv = solve!(L, Ladj, b, g, x, args...)
	return x, slv
end

# when linear operator is composed with g (instead of the least squares term)
# then matrix/linear operator arguments is passed after g in the arguments list
# in this case the dual problem is solved, then the primal solution is recovered

function solve(b::AbstractArray, g::ProximableFunction, A, args...)
	y, slv = solve(A', b, Conjugate(g), args...)
	return -A'*y+b, slv
end

function solve(b::AbstractArray, g::ProximableFunction, L::Function, Ladj::Function, args...)
	y, slv = solve(Ladj, L, b, Conjugate(g), args...)
	return -Ladj(y)+b, slv
end

end
