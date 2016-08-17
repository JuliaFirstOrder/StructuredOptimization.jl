VERSION >= v"0.4.0-dev+6521" && __precompile__()

module RegLS

using LinearOperators

RealOrComplexArray = Union{Array{Float64}, Array{Complex{Float64}}}
MatrixLike = Union{AbstractMatrix,LinearOperator}

include("prox.jl")
include("utils.jl")
include("LBFGS.jl")

abstract Solver

include("ista.jl")
include("fista.jl")
include("zerofpr.jl")

export PG,
       FPG,
       ZeroFPR,
       solve,
       normL2,
       normL2sqr,
       normL1,
       normL21,
       normL0,
       indBallL0,
       indBallRank,
       indBox,
       indBallInf

solve(A::MatrixLike, args...) = solve(x -> A*x, y -> A'*y, args...)

solve(L::Function, Ladj::Function, b::Array, g::Function, x::Array) =
  solve(L, Ladj, b, g, x, ZeroFPR())

end
