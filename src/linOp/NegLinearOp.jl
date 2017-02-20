import Base: +, -

immutable NegLinearOp{D1,D2} <: LinearOp{D1,D2}
	x::OptVar
	A::LinearOp{D1,D2}
end
size(A::NegLinearOp) = size(A.A)

+(A::LinearOp) =  A
-{D1,D2}(A::LinearOp{D1,D2}) =  NegLinearOp{D1,D2}(A.x,A)
-(S::NegLinearOp) =  S.A
+(S::NegLinearOp) =  S

fun_name(S::NegLinearOp) =  "- "*fun_name(S.A)

transpose(S::NegLinearOp) = NegLinearOp(S.x, S.A')

function A_mul_B!(y::AbstractArray,S::NegLinearOp,b::AbstractArray)  
	A_mul_B!(y, S.A, b)
	y .= -y 
end
