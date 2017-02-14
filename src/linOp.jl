abstract LinearOp{D1<: Union{Real,Complex},D2<: Union{Real,Complex}}
#abstract DiagGramOp <: LinearOp
#abstract DiagonalOp <: DiagGramOp


import Base: *,  
             A_mul_B!, 
						 Ac_mul_B!, 
						 transpose,
						 size

immutable OptVar{T<:Union{Real,Complex}}
	x::AbstractArray{T}
end
size(x::OptVar) = size(x.x)
size(A::LinearOp) = (size(A.x),size(A.y))

export OptVar

Ac_mul_B!(y::AbstractArray,A::LinearOp,b::AbstractArray)  = A_mul_B!(y, A',b)


include("linOp/FullOp.jl")
include("linOp/Reshape.jl")
include("linOp/NestedLinearOp.jl")
include("linOp/DFT.jl")
include("linOp/DCT.jl")


function Base.show(io::IO, f::LinearOp)
  println(io, "description : ", fun_name(f))
  println(io, "domain      : ", fun_dom(f))
end

fun_name(  f) = "n/a"
fun_dom(   f) = "n/a"

fun_dom{D1<:Complex, D2<:Real   }(A::LinearOp{D1,D2}) = "ℂ^$(size(A)[1]) →  ℝ^$(size(A)[2])"
fun_dom{D1<:Complex, D2<:Complex}(A::LinearOp{D1,D2}) = "ℂ^$(size(A)[1]) →  ℂ^$(size(A)[2])"
fun_dom{D1<:Real,    D2<:Complex}(A::LinearOp{D1,D2}) = "ℝ^$(size(A)[1]) →  ℂ^$(size(A)[2])"
fun_dom{D1<:Real,    D2<:Real   }(A::LinearOp{D1,D2}) = "ℝ^$(size(A)[1]) →  ℝ^$(size(A)[2])"

fun_dom{D1<:Complex}(A::LinearOp{D1}) = "ℂ^$(size(A)[1]) →  ℂ^$(size(A)[2])"
fun_dom{D1<:Real   }(A::LinearOp{D1}) = "ℝ^$(size(A)[1]) →  ℝ^$(size(A)[2])"
