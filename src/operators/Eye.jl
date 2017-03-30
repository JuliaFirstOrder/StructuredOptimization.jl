import Base: eye
export Eye

immutable Eye{D1,D2} <: IdentityOperator{D1,D2}
	sign::Bool
	dim::Tuple

	Eye(sign,dim) = new(sign,dim)
	Eye(dim) = new(true,dim)
end
size(A::Eye) = (A.dim ,A.dim )
-{D1,D2}(A::Eye{D1,D2}) = Eye{D1,D2}(false == sign(A), A.dim) 

Eye(T::Type, dim::Tuple        ) = Eye{T,T}(dim)
Eye(T::Type, dim::Vararg{Int64}) = Eye{T,T}(dim)
Eye(dim::Vararg{Int64}) = Eye(Float64, dim)
Eye(dim::Tuple        ) = Eye(Float64, dim)

eye{D1}(x::Variable{D1}) = Affine([x], Eye(D1,size(x)), Eye(D1,size(x)),
				Nullable{AbstractArray}() )

function uA_mul_B!{T}(y::AbstractArray{T},A::Eye,b::AbstractArray{T})
	y .= b 
end

transpose{D1}(A::Eye{D1,D1} ) = A
inv{D1}(A::Eye{D1,D1} )       = A

fun_name(A::Eye)  = "Identity Operator"

eye(B::AffineOperator, args...) = NestedLinearOperator(eye, B, args...)

isInvertable(A::Eye) = true
