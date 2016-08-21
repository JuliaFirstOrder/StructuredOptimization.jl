# ------------------------------------------------------------------------------
# prox.jl - library of nonsmooth functions and associated proximal mappings
# ------------------------------------------------------------------------------

abstract ProximableFunction

include("prox/normL2.jl")
include("prox/normL2sqr.jl")
include("prox/normL1.jl")
include("prox/normL21.jl")
include("prox/normL0.jl")
include("prox/elasticNet.jl")
include("prox/indAffine.jl")
include("prox/indBallL0.jl")
include("prox/indBallL2.jl")
include("prox/indBallL20.jl")
include("prox/indBallRank.jl")
include("prox/indBox.jl")

function Base.show(io::IO, f::ProximableFunction)
  println("type       : ", fun_type(f))
  println("expression : ", fun_expr(f))
  print(  "parameters : ", fun_params(f))
end

fun_type(f::ProximableFunction) = "n/a"
fun_expr(f::ProximableFunction) = "n/a"
fun_params(f::ProximableFunction) = "n/a"

return
