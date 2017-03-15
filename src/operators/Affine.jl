import Base: +, -

immutable Affine{Op<:LinearOperator } <: AffineOperator
	A::Op
	b::AbstractArray
end
variable(A::Affine) = variable(A.A)

fun_name(A::Affine) = "Affine "*fun_name(A.A)
fun_dom(A::Affine)  = fun_dom(A.A)

+(A::LinearOperator,b::AbstractArray) = Affine(A, b)
-(A::LinearOperator,b::AbstractArray) = Affine(A,-b)
+(b::AbstractArray, A::LinearOperator) = Affine(A, b)
-(b::AbstractArray, A::LinearOperator) = Affine(-A, b)

+(x::OptVar,b::AbstractArray) = Affine(eye(x), b)
-(x::OptVar,b::AbstractArray) = Affine(eye(x),-b)
+(b::AbstractArray,x::OptVar) = Affine(eye(x), b)
-(b::AbstractArray,x::OptVar) = Affine(-eye(x), b)

function A_mul_B!(y::AbstractArray, A::Affine, x::AbstractArray) 
	A_mul_B!(y,A.A,x)
	y .+= A.b 
end

function *(A::Affine,x::AbstractArray) 
	y1  = A.A*x
	y1 .+= A.b 
end

transpose(A::Affine) = A.A'

