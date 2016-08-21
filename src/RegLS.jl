VERSION >= v"0.4.0-dev+6521" && __precompile__()

module RegLS

using LinearOperators

RealOrComplexArray = Union{Array{Float64}, Array{Complex{Float64}}}
MatrixLike = Union{AbstractMatrix,LinearOperator}

include("prox.jl")
include("utils.jl")
include("lbfgs.jl")

abstract Solver
abstract ForwardBackwardSolver <: Solver

include("pg.jl")
include("fpg.jl")
include("zerofpr.jl")

export PG,
       FPG,
       ZeroFPR,
       solve,
       prox,
       normL2,
       normL2sqr,
       normL1,
       normL21,
       normL0,
       indBallL0,
       indBallL2,
       indBallRank,
       indBox,
       indBallInf,
       indBallL20,
       elasticNet,
       indAffine

solve(A::MatrixLike, b::Array, g::ProximableFunction) =
  solve(x -> A*x, y -> A'*y, b, g, zeros(size(A,2)))

solve(A::MatrixLike, b::Array, g::ProximableFunction, args...) =
  solve(x -> A*x, y -> A'*y, b, g, args...)

solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array) =
  solve(L, Ladj, b, g, x, ZeroFPR())

end
