# ------------------------------------------------------------------------------
# prox.jl - library of nonsmooth functions and associated proximal mappings
# ------------------------------------------------------------------------------

abstract ProximableFunction

include("prox/normL2.jl")
include("prox/normL2sqr.jl")
include("prox/normL1.jl")
include("prox/normL21.jl")
include("prox/normL0.jl")
include("prox/indBallL0.jl")
include("prox/indBallL20.jl")
include("prox/indBallRank.jl")
include("prox/indBox.jl")

return
