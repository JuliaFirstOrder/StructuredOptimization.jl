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

end
