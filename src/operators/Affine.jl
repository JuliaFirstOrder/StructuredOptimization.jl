import Base: +, -

immutable Affine{Op<:LinearOperator, V <:Union{OptVar, Array{OptVar,1}} } <: AffineOperator
	x::V
	A::Op
	b::AbstractArray
end

fun_name(A::Affine) = "Affine "*fun_name(A.A)
fun_dom(A::Affine)  = fun_dom(A.A)

+(A::LinearOperator,b::AbstractArray) = Affine(A.x, A, b)
-(A::LinearOperator,b::AbstractArray) = Affine(A.x, A,-b)

+(x::OptVar,b::AbstractArray) = Affine(x,eye(x), b)
-(x::OptVar,b::AbstractArray) = Affine(x,eye(x),-b)

function A_mul_B!(y::AbstractArray, A::Affine, x::AbstractArray) 
	A_mul_B!(y,A.A,x)
	y .+= A.b 
end

function *(A::Affine,x::AbstractArray) 
	y1  = A.A*x
	y1 .+= A.b 
end

transpose(A::Affine) = A.A'

