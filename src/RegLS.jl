__precompile__()

module RegLS

using Prox

RealOrComplexArray = Union{Array{Float64}, Array{Complex{Float64}}}
RealOrComplexVector = Union{Array{Float64,1}, Array{Complex{Float64},1}}
RealOrComplexMatrix = Union{Array{Float64,2}, Array{Complex{Float64},2}}

include("solvers.jl")

export solve

solve(A, b::Array, g::ProximableFunction) =
  solve(x -> A*x, y -> A'*y, b, g, zeros(size(A,2)))

solve(A, b::Array, g::ProximableFunction, args...) =
  solve(x -> A*x, y -> A'*y, b, g, args...)

solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array) =
  solve(L, Ladj, b, g, x, ZeroFPR())

end
