abstract LinearOp
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
export OptVar

Ac_mul_B!(y::AbstractArray,A::LinearOp,b::AbstractArray)  = A_mul_B!(y, A',b)


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

function fun_dom(A::LinearOp)
	a = A.dom[1] <: Complex ? "ℂ^$(A.dim[1])   " : "R^$(A.dim[1])   "
	b = A.dom[2] <: Complex ? "→  ℂ^$(A.dim[2])" : "→  R^$(A.dim[2])"
	return a*b
end
