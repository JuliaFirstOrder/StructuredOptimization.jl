abstract AffineOperator
abstract LinearOperator{D1,D2} <: AffineOperator
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
  ==

size(A::LinearOperator, i::Int64) = size(A)[i]
ndims(A::LinearOperator) = (length(size(A,1)),length(size(A,2)))
ndims(A::LinearOperator, i::Int64) = length(size(A,i))

Ac_mul_B!(y::AbstractArray,A::LinearOperator,b::AbstractArray)  = A_mul_B!(y, A',b)

function *{D1,D2}(A::LinearOperator{D1,D2},b::AbstractArray)
	y = Array{D2}(size(A,2))
	A_mul_B!(y,A,b)
	return y
end

.*(A::LinearOperator,b) = A*b

include("operators/utils.jl")
include("operators/Affine.jl")
include("operators/FullOp.jl")
include("operators/Reshape.jl")
include("operators/NestedLinearOp.jl")
include("operators/DFT.jl")
include("operators/DCT.jl")
include("operators/DiagOp.jl")
include("operators/Eye.jl")
include("operators/GetIndex.jl")
include("operators/Zeros.jl")
include("operators/HCAT.jl")
include("operators/VCAT.jl")
include("operators/Plus.jl")
include("operators/LBFGS.jl")

function Base.show{Op <: AffineOperator }(io::IO, f::Op)
  println(io, "description : ", fun_name(f))
  println(io, "domain      : ", fun_dom(f))
end

fun_name(  f) = "n/a"
fun_dom(   f) = "n/a"
fun_par(   f) = "n/a"

fun_dom{D1<:Complex, D2<:Real   }(A::LinearOperator{D1,D2}) = "ℂ^$(size(A,1)) →  ℝ^$(size(A,2))"
fun_dom{D1<:Complex, D2<:Complex}(A::LinearOperator{D1,D2}) = "ℂ^$(size(A,1)) →  ℂ^$(size(A,2))"
fun_dom{D1<:Real,    D2<:Complex}(A::LinearOperator{D1,D2}) = "ℝ^$(size(A,1)) →  ℂ^$(size(A,2))"
fun_dom{D1<:Real,    D2<:Real   }(A::LinearOperator{D1,D2}) = "ℝ^$(size(A,1)) →  ℝ^$(size(A,2))"

isEye(A::AffineOperator) = typeof(A.A) <: IdentityOperator 
isEye(A::LinearOperator) = typeof(A) <: IdentityOperator 

isDiagonal(A::AffineOperator) = typeof(A.A) <: DiagonalOperator 
isDiagonal(A::LinearOperator) = typeof(A) <: DiagonalOperator 

isAbsorbable(A::AffineOperator) = typeof(A.A) <: DiagonalOperator 
isAbsorbable(A::LinearOperator) = typeof(A) <: DiagonalOperator 
#this will be changed with typeof(A) <: GramDiagonal 

isInvertable(A::AffineOperator) = isInvertable(A.A) 
isInvertable(A::LinearOperator) = false

inv(A::Affine) = inv(A.A)

==(A::Affine        , B::Affine        ) = A.A == B.A 
==(A::Affine        , B::LinearOperator) = A.A ==   B 
==(A::LinearOperator, B::Affine        ) = A   == B.A 







