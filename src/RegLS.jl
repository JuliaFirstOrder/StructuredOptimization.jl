VERSION >= v"0.4.0-dev+6521" && __precompile__()

module RegLS

RealOrComplexArray = Union{Array{Float64}, Array{Complex{Float64}}}

include("prox.jl")
include("solvers.jl")

export solve

solve(A, b::Array, g::ProximableFunction) =
  solve(x -> A*x, y -> A'*y, b, g, zeros(size(A,2)))

solve(A, b::Array, g::ProximableFunction, args...) =
  solve(x -> A*x, y -> A'*y, b, g, args...)

solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array) =
  solve(L, Ladj, b, g, x, ZeroFPR())

end
