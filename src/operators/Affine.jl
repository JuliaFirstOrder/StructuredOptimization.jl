import Base: +, -

immutable Affine{Op<:LinearOp, V <:Union{OptVar, Array{OptVar,1}} } <: AffineOp 
	x::V
	A::Op
	b::AbstractArray
	pos::Bool
end

fun_name(A::Affine) = "Affine "*fun_name(A.A)
fun_dom(A::Affine)  = fun_dom(A.A)

+(A::LinearOp,b::AbstractArray) = Affine(A.x, A, b, true )
-(A::LinearOp,b::AbstractArray) = Affine(A.x, A, b, false)

function A_mul_B!(y::AbstractArray, A::Affine, x::AbstractArray) 
	A_mul_B!(y,A.A,x)
	A.pos ? y .+= A.b : y .-= A.b
end

function *(A::Affine,x::AbstractArray) 
	y1  = A.A*x
	A.pos ? y1 .+= A.b : y1 .-=  A.b
end

transpose(A::Affine) = A.A'

