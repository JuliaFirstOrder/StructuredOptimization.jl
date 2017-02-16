import Base: eye

immutable Eye{D1,D2} <: LinearOp{D1,D2}
	dim::Tuple
end
size(A::Eye) = A.dim

eye{D1}(x::OptVar{D1}) = Eye{D1,D1}((size(x),size(x)))
eye{D1}(x::OptVar{D1},n::Int64) = Eye{D1,D1}((size(x),(n,)))

function *{T}(A::Eye, b::AbstractArray{T})  
	y = zeros(T,size(A,2)[1])
	A_mul_B!(y,A,b)
	return y
end

function A_mul_B!{T}(y::AbstractArray{T},A::Eye,b::AbstractArray{T})
	n = min(length(b), size(A,2)[1])
	y .= (*).(y,0)
	copy!(y,1,b,1,n)
end

transpose{D1}(A::Eye{D1,D1} ) = Eye{D1,D1}((A.dim[2],A.dim[1]))

fun_name(A::Eye)  = "Identity Operator"

eye(B::LinearOp, args...) = NestedLinearOp(eye,B, args...)
