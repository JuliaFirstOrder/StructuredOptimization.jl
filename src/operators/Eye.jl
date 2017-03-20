import Base: eye
export Eye

immutable Eye{D1,D2} <: IdentityOperator{D1,D2}
	dim::Tuple
end
size(A::Eye) = (A.dim ,A.dim )

Eye(T::Type, dim...) = Eye{T,T}(dim...)
Eye(dim...) = Eye(Float64, dim)

eye{D1}(x::OptVar{D1}) = Affine([x], Eye(D1,size(x)), Eye(D1,size(x)),
				Nullable{Vector{AbstractArray}}() )

function A_mul_B!{T}(y::AbstractArray{T},A::Eye,b::AbstractArray{T})
	y .= b
end

transpose{D1}(A::Eye{D1,D1} ) = A
inv{D1}(A::Eye{D1,D1} )       = A

fun_name(A::Eye)  = "Identity Operator"

eye(B::AffineOperator, args...) = NestedLinearOperator(eye, B, args...)

isInvertable(A::Eye) = true
