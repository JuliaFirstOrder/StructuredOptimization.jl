struct Expression{N,A<:AbstractOperator} <: AbstractExpression
  x::NTuple{N,Variable}
  L::A
  function Expression{N}(x::NTuple{N,Variable}, L::A) where {N,A<:AbstractOperator}
    # checks on L
    ndoms(L,1) > 1 && throw(ArgumentError(
      "Cannot create expression with LinearOperator with `ndoms(L,1) > 1`"
     ))
    #checks on x
    szL = size(L,2)
    szx = size.(x)
    check_sz = length(szx) == 1 ? szx[1] != szL : szx != szL
    check_sz && throw(ArgumentError(
      "Size of the operator domain $(size(L, 2)) must match size of the variable $(size.(x))"
     ))
    dmL = domainType(L)
    dmx = eltype.(x)
    check_dm = length(dmx) == 1 ? dmx[1] != dmL : dmx != dmL
    check_dm && throw(ArgumentError(
      "Type of the operator domain $(domainType(L)) must match type of the variable $(eltype.(x))"
     ))
    new{N,A}(x,L)
  end
end

struct AdjointExpression{E <: AbstractExpression} <: AbstractExpression
  ex::E
end

import Base: adjoint

adjoint(ex::AbstractExpression) = AdjointExpression(convert(Expression,ex))
adjoint(ex::AdjointExpression) = ex.ex

include("utils.jl")
include("multiplication.jl")
include("addition.jl")
include("abstractOperator_bind.jl")
