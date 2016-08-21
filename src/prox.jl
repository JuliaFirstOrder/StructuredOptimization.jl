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
include("prox/indSOC.jl")

function call(f::ProximableFunction, x::Array)
  error("call is not implemented for type ", typeof(f))
end

"""
  prox(f::ProximableFunction, γ::Float64, x::Array)

Computes the proximal point of `x` with respect to function `f`
and parameter `γ > 0`, that is

  y = argmin_z { f(z) + 1/(2γ)||z-x||^2 }

and returns `y` and `f(y)`.
"""

function prox(f::ProximableFunction, gamma::Float64, x::Array)
  error("prox is not implemented for type ", typeof(f))
end

function Base.show(io::IO, f::ProximableFunction)
  println(io, "name       : ", fun_name(f))
  println(io, "type       : ", fun_type(f))
  println(io, "expression : ", fun_expr(f))
  print(  io, "parameters : ", fun_params(f))
end

fun_name(f::ProximableFunction) = "n/a"
fun_type(f::ProximableFunction) = "n/a"
fun_expr(f::ProximableFunction) = "n/a"
fun_params(f::ProximableFunction) = "n/a"

return
