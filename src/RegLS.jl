__precompile__()

module RegLS

using ProximalOperators

include("solvers.jl")

export solve, solve!

function solve(A, b::AbstractArray, g::ProximableFunction, args...)

	y = zeros(typeof(b[1]),size(A,1))
	y2 = zeros(typeof(b[1]),size(A,2))

	L! = x -> A_mul_B!(y, A, x)
	Ladj! = x -> Ac_mul_B!(y2, A, x)

	solve(L!, Ladj!, b, g, args...)

end

solve(A, b::AbstractArray, g::ProximableFunction) =
	solve(A, b::AbstractArray, g::ProximableFunction, zeros(eltype(b), size(A,2)))

solve(L::Function, Ladj::Function, b::AbstractArray, g::ProximableFunction, x::Array) =
  solve(L, Ladj, b, g, x, ZeroFPR())

# when linear operator is composed with g (instead of the least squares term)

solve(b::AbstractArray, g::ProximableFunction, A, args...) =
	solve(A', b, Conjugate(g), args...)
solve(b::AbstractArray, g::ProximableFunction, L::Function, Ladj::Function, x::Array) =
	solve(Ladj, L, b, Conjugate(g), x, ZeroFPR())

end
