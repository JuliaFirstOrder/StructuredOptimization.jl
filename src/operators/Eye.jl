export Eye

immutable Eye <: IdentityOperator
	domainType::Type
	dim::Tuple
end

size(L::Eye) = (L.dim, L.dim)
  domainType(L::Eye) = L.domainType
codomainType(L::Eye) = L.domainType

Eye{T}(x::AbstractArray{T}) = Eye(T,size(x))

A_mul_B!{T}(y::AbstractArray{T}, L::Eye, b::AbstractArray{T}) = y .= b
Ac_mul_B!(y, L::Eye, b) = A_mul_B!(y,L,b) 

transpose(L::Eye) = L
inv(L::Eye ) = L

fun_name(L::Eye)  = "Identity Operator"
isInvertible(L::Eye) = true
