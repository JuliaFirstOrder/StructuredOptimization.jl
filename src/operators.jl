abstract AffineOperator
abstract LinearOperator{D1,D2}
abstract DiagonalOperator{D1,D2} <: LinearOperator{D1,D2}
abstract IdentityOperator{D1,D2} <: DiagonalOperator{D1,D2}
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

size(A::LinearOperator, i) = size(A)[i]
nblocks(A::LinearOperator) = (length(size(A,1)),length(size(A,2)))
nblocks(A::LinearOperator, i) = length(size(A,i))

#   domainType{D1,D2}(A::LinearOperator{D1,D2}) = D1
# codomainType{D1,D2}(A::LinearOperator{D1,D2}) = D2

# Ac_mul_B!(y::AbstractArray,A::LinearOperator,b::AbstractArray)  = A_mul_B!(y, A',b)
#  A_mul_B!(y::AbstractArray,A::LinearOperator,b::AbstractArray) =
# sign(A) ? uA_mul_B!(y,A,b) : uA_mul_B!(y,A,-b)

# function *{D1,D2}(A::LinearOperator{D1,D2},b::AbstractArray)
# 	y = zeros(D2,size(A,2))
# 	A_mul_B!(y,A,b)
# 	return y
# end

# .*(A::LinearOperator,b) = A*b
#
# +(A::LinearOperator) = A
# sign(A::LinearOperator) = A.sign ? true : false

function (*)(L::LinearOperator, x::AbstractArray)
  if nblocks(L, 1) == 1
    y = zeros(size(L, 1))
  else
    y = Vector()
    for s in size(L, 1)
      push!(y, zeros(s))
    end
  end
  A_mul_B!(y, L, x)
end

include("operators/Affine.jl")
include("operators/Eye.jl")
include("operators/MatrixOp.jl")
include("operators/Reshape.jl")
include("operators/Compose.jl")
include("operators/DFT.jl")
include("operators/FiniteDiff.jl")
include("operators/TV.jl")
include("operators/DCT.jl")
include("operators/Conv.jl")
include("operators/DiagOp.jl")
include("operators/GetIndex.jl")
include("operators/Empty.jl")
include("operators/HCAT.jl")
include("operators/VCAT.jl")
include("operators/Sum.jl")
include("operators/Scale.jl")
include("operators/Transpose.jl")
include("operators/LBFGS.jl")
include("operators/Zeros.jl")
include("operators/utils.jl")

function Base.show(io::IO, L::LinearOperator)
  println(io, "description : ", fun_name(L))
  print(  io, "type        : ", fun_type(L))
end

fun_name(  L) = "n/a"
fun_type(  L) = "$(size(L,2)) → $(size(L,1))"
fun_par(   L) = "n/a"

# fun_domain{D1<:Complex, D2}(A::LinearOperator{D1,D2})   = "ℂ^$(size(A,2))"
# fun_domain{D1<:Real,    D2}(A::LinearOperator{D1,D2})   = "ℝ^$(size(A,2))"
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
