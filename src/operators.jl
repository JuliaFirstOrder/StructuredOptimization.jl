abstract AffineOperator
abstract LinearOperator
abstract DiagonalOperator <:   LinearOperator
abstract IdentityOperator <: DiagonalOperator
#abstract DiagonalOp <: DiagGramOp

import Base:
  *,
  .*,
  A_mul_B!,
  Ac_mul_B!,
  transpose,
  inv,
  size,
  ndims,
  +,
  -,
  ==


+(L::LinearOperator) = L

function *(L::LinearOperator,b::AbstractArray)
	y = zeros(codomainType(L),size(L,1))        
	A_mul_B!(y,L,b)        
	return y 
end

#calcolus
include("operators/Scale.jl")
include("operators/Transpose.jl")
include("operators/Sum.jl")
include("operators/Reshape.jl")
include("operators/Compose.jl")
include("operators/HCAT.jl")
include("operators/VCAT.jl")

include("operators/Zeros.jl")
include("operators/Eye.jl")
include("operators/DiagOp.jl")
include("operators/GetIndex.jl")
include("operators/MatrixOp.jl")
include("operators/DFT.jl")
include("operators/DCT.jl")
include("operators/FiniteDiff.jl")
include("operators/Variation.jl")
include("operators/Conv.jl")
include("operators/Xcorr.jl")
include("operators/LBFGS.jl")
include("operators/utils.jl")

size(L::LinearOperator, i) = size(L)[i]

#usually domain is preserved, if not one has to redefine these functions
  domainType(L::LinearOperator) = L.domainType
codomainType(L::LinearOperator) = L.domainType

isEye(L::LinearOperator) = typeof(L) <: IdentityOperator
isDiagonal(L::LinearOperator) = typeof(L) <: DiagonalOperator
# this will be changed with typeof(L) <: GramDiagonal
isGramDiagonal(L::LinearOperator) = typeof(L) <: DiagonalOperator

isInvertible(L::LinearOperator) = false

function Base.show(io::IO, L::LinearOperator)
  println(io, "description : ", fun_name(L))
  print(  io, "type        : ", fun_type(L))
end

fun_type(  L) = fun_domain(L)*" → "*fun_codomain(L)

fun_domain(L::LinearOperator)   =   domainType(L) <: Complex ? "ℂ^$(size(L,2))" : "ℝ^$(size(L,2))"
fun_codomain(L::LinearOperator) = codomainType(L) <: Complex ? "ℂ^$(size(L,1))" : "ℝ^$(size(L,1))"
