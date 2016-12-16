__precompile__()

module RegLS

using ProximalOperators

include("solvers.jl")

export solve, solve!

function solve(A, b::Array, g::ProximableFunction, args...)

	y = zeros(typeof(b[1]),size(A,1))
	y2 = zeros(typeof(b[1]),size(A,2))

	L! = x -> A_mul_B!(y, A, x)
	Ladj! = x -> Ac_mul_B!(y2, A, x)

	solve(L!, Ladj!, b, g, args...)

end

function solve(A, b::Array, g::ProximableFunction)

	y = zeros(typeof(b[1]),size(A,1))
	y2 = zeros(typeof(b[1]),size(A,2))

	L! = x -> A_mul_B!(y, A, x)
	Ladj! = x -> Ac_mul_B!(y2, A, x)

	solve(L!, Ladj!, b, g, zeros(typeof(b[1]),size(A,2)))

end

solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array) =
  solve(L, Ladj, b, g, x, ZeroFPR())

# when linear operator is composed with g (instead of the least squares term)

solve(b::Array, g::ProximableFunction, A) = solve(A', b, Conjugate(g))
solve(b::Array, g::ProximableFunction, A, args...) = solve(A', b, Conjugate(g), args...)
solve(b::Array, g::ProximableFunction, L::Function, Ladj::Function, x::Array) =
  solve(Ladj, L, b, Conjugate(g), x, ZeroFPR())

end
