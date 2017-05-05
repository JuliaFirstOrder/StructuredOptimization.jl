export Conv

immutable Conv{N} <: LinearOperator
	domainType::Type
	dim_in::NTuple{N,Int}
	h::AbstractVector
end

size(L::Conv) = (L.dim_in[1]+length(L.h)-1,), L.dim_in

# Constructors

Conv{D1}(dim_in::Int,  h::AbstractVector{D1}) = Conv(eltype(h),(dim_in,), h)
Conv{D1}(x::AbstractVector{D1}, h::AbstractVector{D1}) = Conv(eltype(x),size(x), h)

# Operators

function A_mul_B!{T}(y::AbstractVector{T},A::Conv,b::AbstractVector{T})
		y .= conv(A.h,b)
end

function Ac_mul_B!{T}(y::AbstractVector{T},A::Conv,b::AbstractVector{T})
		y .= xcorr(b,A.h)[size(A,1)[1]:end-length(A.h)+1]
end

# Properties
fun_name(A::Conv)  = "Convolution Operator"

