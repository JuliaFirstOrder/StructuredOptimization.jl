import Base: eye

immutable Eye{D1,D2} <: IdentityOperator{D1,D2}
	x::OptVar
end
size(A::Eye) = (size(A.x),size(A.x))

eye{D1}(x::OptVar{D1}) = Eye{D1,D1}(x)

function A_mul_B!{T}(y::AbstractArray{T},A::Eye,b::AbstractArray{T})
	y .= b
end

transpose{D1}(A::Eye{D1,D1} ) = A
inv{D1}(A::Eye{D1,D1} )       = A

fun_name(A::Eye)  = "Identity Operator"

eye(B::LinearOperator, args...) = NestedLinearOperator(eye, B, args...)

isInvertable(A::Eye) = true
