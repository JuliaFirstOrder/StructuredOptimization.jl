abstract type LinearOperator end

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


include("operators/MyOperator.jl")
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
include("operators/Filt.jl")
include("operators/MIMOFilt.jl")
include("operators/ZeroPad.jl")
include("operators/Xcorr.jl")
include("operators/LBFGS.jl")
include("operators/BlkDiagLBFGS.jl")
include("operators/utils.jl")

#calcolus
include("operators/Scale.jl")
include("operators/Transpose.jl")
include("operators/Sum.jl")
include("operators/Reshape.jl")
include("operators/Compose.jl")
include("operators/HCAT.jl")
include("operators/VCAT.jl")

size(L::LinearOperator, i::Int) = size(L)[i]
ndims(L::LinearOperator)   = length(size(L,1)),length(size(L,2))
ndims(L::LinearOperator, i::Int)   = ndims(L)[i]

#usually domain is preserved, if not one has to redefine these functions
  domainType(L::LinearOperator) = L.domainType
codomainType(L::LinearOperator) = L.domainType

isScaled(L::LinearOperator)         = false

isEye(L::LinearOperator)            = false
isDiagonal(L::LinearOperator)       = false
isGramDiagonal(L::LinearOperator)   = false

isInvertible(L::LinearOperator)     = false
isFullRowRank(L::LinearOperator)    = false
isFullColumnRank(L::LinearOperator) = false

function Base.show(io::IO, L::LinearOperator)
  println(io, "description : ", fun_name(L))
  print(  io, "type        : ", fun_type(L))
end

fun_type(  L) = fun_domain(L)*" → "*fun_codomain(L)

fun_domain(L::LinearOperator)   =   domainType(L) <: Complex ? "ℂ^$(size(L,2))" : "ℝ^$(size(L,2))"
fun_codomain(L::LinearOperator) = codomainType(L) <: Complex ? "ℂ^$(size(L,1))" : "ℝ^$(size(L,1))"
