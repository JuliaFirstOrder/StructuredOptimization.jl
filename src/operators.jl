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
  sign,
  ==

size(L::LinearOperator, i) = size(L)[i]
nblocks(L::LinearOperator) = (length(size(L,1)),length(size(L,2)))
nblocks(L::LinearOperator, i) = length(size(L,i))

  domainType(L::LinearOperator) =   L.domainType
codomainType(L::LinearOperator) = L.codomainType


# function *{D1,D2}(A::LinearOperator{D1,D2},b::AbstractArray)
# 	y = zeros(D2,size(A,2))
# 	A_mul_B!(y,A,b)
# 	return y
# end

+(L::LinearOperator) = L

function *(L::LinearOperator,b::AbstractArray)
	y = zeros(codomainType(L),size(L,1))        
	A_mul_B!(y,L,b)        
	return y 
end





include("operators/Eye.jl")
#include("operators/MatrixOp.jl")
#include("operators/Reshape.jl")
#include("operators/Compose.jl")
#include("operators/DFT.jl")
#include("operators/FiniteDiff.jl")
#include("operators/TV.jl")
#include("operators/DCT.jl")
#include("operators/Conv.jl")
#include("operators/DiagOp.jl")
#include("operators/GetIndex.jl")
#include("operators/Empty.jl")
#include("operators/HCAT.jl")
#include("operators/VCAT.jl")
#include("operators/Sum.jl")
#include("operators/Scale.jl")
#include("operators/Transpose.jl")
#include("operators/LBFGS.jl")
#include("operators/Zeros.jl")
include("operators/utils.jl")

function Base.show(io::IO, L::LinearOperator)
  println(io, "description : ", fun_name(L))
  print(  io, "type        : ", fun_type(L))
end

fun_name(  L) = "n/a"
fun_type(  L) = fun_domain(L)*" → "*fun_codomain(L)
fun_par(   L) = "n/a"

fun_domain(L::LinearOperator)   =   domainType(L) <: Complex ? "ℂ^$(size(L,2))" : "ℝ^$(size(L,2))"
fun_codomain(L::LinearOperator) = codomainType(L) <: Complex ? "ℂ^$(size(L,1))" : "ℝ^$(size(L,1))"
#
# fun_codomain{D1, D2<:Complex}(A::LinearOperator{D1,D2}) = "ℂ^$(size(A,1))"
# fun_codomain{D1, D2<:Real   }(A::LinearOperator{D1,D2}) = "ℝ^$(size(A,1))"

# fun_type{D1<:Complex, D2<:Real   }(A::LinearOperator{D1,D2}) = "ℂ^$(size(A,1)) →  ℝ^$(size(A,2))"
# fun_type{D1<:Complex, D2<:Complex}(A::LinearOperator{D1,D2}) = "ℂ^$(size(A,1)) →  ℂ^$(size(A,2))"
# fun_type{D1<:Real,    D2<:Complex}(A::LinearOperator{D1,D2}) = "ℝ^$(size(A,1)) →  ℂ^$(size(A,2))"
# fun_type{D1<:Real,    D2<:Real   }(A::LinearOperator{D1,D2}) = "ℝ^$(size(A,1)) →  ℝ^$(size(A,2))"

isEye(A::LinearOperator) = typeof(A) <: IdentityOperator
isDiagonal(A::LinearOperator) = typeof(A) <: DiagonalOperator
# this will be changed with typeof(A) <: GramDiagonal
isGramDiagonal(A::LinearOperator) = typeof(A) <: DiagonalOperator
# is this needed?
isInvertible(A::LinearOperator) = false
