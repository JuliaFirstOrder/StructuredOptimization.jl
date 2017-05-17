export Xcorr

immutable Xcorr{N} <: LinearOperator
	domainType::Type
	dim_in::NTuple{N,Int}
	h::AbstractVector
end

# Constructors

Xcorr{D1}(dim_in::Int,  h::AbstractVector{D1}) = Xcorr(eltype(h),(dim_in,), h)
Xcorr{D1}(x::AbstractVector{D1}, h::AbstractVector{D1}) = Xcorr(eltype(x),size(x), h)

# Mappings

function A_mul_B!{T}(y::AbstractVector{T},A::Xcorr,b::AbstractVector{T})
		y .= xcorr(b,A.h)
end

function Ac_mul_B!{T}(y::AbstractVector{T},A::Xcorr,b::AbstractVector{T})
		l =floor(Int64,size(A,1)[1]/2)
		idx = l+1:l+length(y)
		y .= conv(b,A.h)[idx]
end

# Properties

size(L::Xcorr) = ( 2*max(L.dim_in[1], length(L.h))-1, ), L.dim_in

fun_name(A::Xcorr)  = "Cross Correlation Operator"
