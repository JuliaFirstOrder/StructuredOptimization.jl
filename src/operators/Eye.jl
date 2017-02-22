import Base: eye

immutable Eye{D1,D2} <: LinearOp{D1,D2}
	x::OptVar
	dim::Tuple
end
size(A::Eye) = A.dim

eye{D1}(x::OptVar{D1}) = Eye{D1,D1}(x,(size(x),size(x)))

function A_mul_B!{T}(y::AbstractArray{T},A::Eye,b::AbstractArray{T})
	y .= b
end

transpose{D1}(A::Eye{D1,D1} ) = Eye{D1,D1}(A.x, (A.dim[2],A.dim[1]))

fun_name(A::Eye)  = "Identity Operator"

eye(B::LinearOp, args...) = NestedLinearOp(eye, B, args...)
