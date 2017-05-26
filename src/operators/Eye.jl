export Eye

immutable Eye{T, N} <: LinearOperator
	dim::NTuple{N, Integer}
end

# Constructors

Eye{N}(t::Type, dims::NTuple{N, Integer}) = Eye{t, N}(dims)
Eye(t::Type, dims::Vararg{Integer}) = Eye{t, length(dims)}(dims)

# Mappings

A_mul_B!{T, N}(y::AbstractArray{T, N}, L::Eye{T, N}, b::AbstractArray{T, N}) = y .= b
Ac_mul_B!{T, N}(y::AbstractArray{T, N}, L::Eye{T, N}, b::AbstractArray{T, N}) = A_mul_B!(y, L, b)

# Properties

domainType{T, N}(L::Eye{T, N}) = T
codomainType{T, N}(L::Eye{T, N}) = T

size(L::Eye) = (L.dim, L.dim)

fun_name(L::Eye) = "Identity"

is_eye(L::Eye) = true
is_diagonal(L::Eye) = true
is_gram_diagonal(L::Eye) = true
is_invertible(L::Eye) = true
is_full_row_rank(L::Eye) = true
is_full_column_rank(L::Eye) = true
