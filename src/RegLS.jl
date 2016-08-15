VERSION >= v"0.4.0-dev+6521" && __precompile__()

module RegLS

RealOrComplexArray = Union{Array{Float64}, Array{Complex{Float64}}}

include("prox.jl")
include("utils.jl")
include("LBFGS.jl")
include("ista.jl")
include("fista.jl")
include("zerofpr.jl")

export ista,
       fista,
       zerofpr,
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

solve = zerofpr

end
