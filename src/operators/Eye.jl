export Eye

immutable Eye{N} <: LinearOperator
	domainType::Type
	dim_in::NTuple{N,Int}
end

# Constructors

Eye(dim_in::Vararg{Int64})  = Eye(dim_in)
Eye(domainType::Type, dim_in::Vararg{Int64})  = Eye(domainType,dim_in)
Eye(dim_in::Tuple)          = Eye(Float64,dim_in)
Eye{T}(x::AbstractArray{T}) = Eye(T,size(x))

# Mappings

A_mul_B!{T}(y::AbstractArray{T}, L::Eye, b::AbstractArray{T}) = y .= b
Ac_mul_B!(y, L::Eye, b) = A_mul_B!(y,L,b)

# Transformations
# inv(L::Eye ) = L

#Properties

size(L::Eye) = (L.dim_in, L.dim_in)

fun_name(L::Eye)  = "Identity Operator"

is_diagonal(L::Eye) = true
is_gram_diagonal(L::Eye) = true
is_invertible(L::Eye) = true
is_full_row_rank(L::Eye) = true
is_full_column_rank(L::Eye) = true
