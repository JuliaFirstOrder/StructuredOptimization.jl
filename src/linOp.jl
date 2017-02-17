abstract LinearOp{D1<: Union{Real,Complex},D2<: Union{Real,Complex}}
#abstract DiagGramOp <: LinearOp
#abstract DiagonalOp <: DiagGramOp

import Base: *,  
             A_mul_B!, 
						 Ac_mul_B!, 
						 transpose,
						 size,
						 ndims

immutable OptVar{T<:Union{Real,Complex}}
	x::AbstractArray{T}
end
size(x::OptVar) = size(x.x)
size(A::LinearOp, i::Int64) = size(A)[i]
ndims(A::LinearOp) = (length(size(A,1)),length(size(A,2)))
ndims(A::LinearOp, i::Int64) = length(size(A,i))

export OptVar

Ac_mul_B!(y::AbstractArray,A::LinearOp,b::AbstractArray)  = A_mul_B!(y, A',b)

function *{D1,D2}(A::LinearOp{D1,D2},b::AbstractArray) 
	y = Array{D2}(size(A,2))
	A_mul_B!(y,A,b)
	return y
end


include("linOp/FullOp.jl")
include("linOp/Reshape.jl")
include("linOp/NestedLinearOp.jl")
include("linOp/DFT.jl")
include("linOp/DCT.jl")
include("linOp/Eye.jl")
include("linOp/DiagOp.jl")
include("linOp/GetIndex.jl")
include("linOp/Plus.jl")


function Base.show(io::IO, f::LinearOp)
  println(io, "description : ", fun_name(f))
  println(io, "domain      : ", fun_dom(f))
end

fun_name(  f) = "n/a"
fun_dom(   f) = "n/a"

fun_dom{D1<:Complex, D2<:Real   }(A::LinearOp{D1,D2}) = "ℂ^$(size(A,1)) →  ℝ^$(size(A,2))"
fun_dom{D1<:Complex, D2<:Complex}(A::LinearOp{D1,D2}) = "ℂ^$(size(A,1)) →  ℂ^$(size(A,2))"
fun_dom{D1<:Real,    D2<:Complex}(A::LinearOp{D1,D2}) = "ℝ^$(size(A,1)) →  ℂ^$(size(A,2))"
fun_dom{D1<:Real,    D2<:Real   }(A::LinearOp{D1,D2}) = "ℝ^$(size(A,1)) →  ℝ^$(size(A,2))"
