export Eye

immutable Eye <: IdentityOperator
	domainType::Type
	dim_in::Tuple
end

size(L::Eye) = (L.dim_in, L.dim_in)

Eye(dim_in::Tuple)             = Eye(Float64,dim_in)
Eye{T}(x::AbstractArray{T}) = Eye(T,size(x))

A_mul_B!{T}(y::AbstractArray{T}, L::Eye, b::AbstractArray{T}) = y .= b
Ac_mul_B!(y, L::Eye, b) = A_mul_B!(y,L,b) 

transpose(L::Eye) = L
inv(L::Eye ) = L

fun_name(L::Eye)  = "Identity Operator"
isInvertible(L::Eye) = true
