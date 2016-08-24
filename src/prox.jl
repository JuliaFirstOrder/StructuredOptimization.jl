# ------------------------------------------------------------------------------
# prox.jl - library of nonsmooth functions and associated proximal mappings
# ------------------------------------------------------------------------------

export prox

export DistL2,
       SqrDistL2,
       ElasticNet,
       IndAffine,
       IndBallInf,
       IndBallL0,
       IndBallL2,
       IndBallL20,
       IndBallRank,
       IndBox,
       IndHalfspace,
       IndNonnegative,
       IndSimplex,
       IndSOC,
       NormL0,
       NormL1,
       NormL2,
       NormL21,
       SqrNormL2

# This hierarchy of abstract types (or a similar one) may be useful.
# Unfortunately Julia does not allow for multiple inheritance.
#
# ProximableFunction --+-- NormFunction
#                      |
#                      +-- IndicatorFunction -- IndicatorConvex
#

abstract ProximableFunction
abstract NormFunction <: ProximableFunction
abstract IndicatorFunction <: ProximableFunction
abstract IndicatorConvex <: IndicatorFunction

include("prox/distL2.jl")
include("prox/elasticNet.jl")
include("prox/normL2.jl")
include("prox/normL1.jl")
include("prox/normL21.jl")
include("prox/normL0.jl")
include("prox/indAffine.jl")
include("prox/indBallL0.jl")
include("prox/indBallL2.jl")
include("prox/indBallL20.jl")
include("prox/indBallRank.jl")
include("prox/indBox.jl")
include("prox/indSOC.jl")
include("prox/indHalfspace.jl")
include("prox/indSimplex.jl")
include("prox/sqrDistL2.jl")
include("prox/sqrNormL2.jl")

function call(f::ProximableFunction, x)
  error("call is not implemented for type ", typeof(f))
end

"""
  prox(f::ProximableFunction, γ::Float64, x::Array)

Computes the proximal point of `x` with respect to function `f`
and parameter `γ > 0`, that is

  y = argmin_z { f(z) + 1/(2γ)||z-x||^2 }

and returns `y` and `f(y)`.
"""

function prox(f::ProximableFunction, gamma::Float64, x)
  error("prox is not implemented for type ", typeof(f))
end

function Base.show(io::IO, f::ProximableFunction)
  println(io, "description : ", fun_name(f))
  println(io, "type        : ", fun_type(f))
  println(io, "expression  : ", fun_expr(f))
  print(  io, "parameters  : ", fun_params(f))
end

fun_name(  f::ProximableFunction) = "n/a"
fun_type(  f::ProximableFunction) = "n/a"
fun_expr(  f::ProximableFunction) = "n/a"
fun_params(f::ProximableFunction) = "n/a"

return
