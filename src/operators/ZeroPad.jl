export ZeroPad

immutable ZeroPad{N} <: LinearOperator
	domainType::Type
	dim_out::NTuple{N,Int}
	dim_in::NTuple{N,Int}
	function ZeroPad{N}(domainType, dim_in, zp::Vararg{Int}) where {N}
		any([zp...].<0) && throw(ArgumentError("zero padding cannot be negative"))
		dim_out = tuple(
		  [i<=length(zp) ? dim_in[i]+zp[i] : dim_in[i] for i in eachindex(dim_in)]...)

		new(domainType,dim_out,dim_in)
	end
end

size(L::ZeroPad) = L.dim_out, L.dim_in

# Constructors
ZeroPad{N}(dim_in::Tuple, zp::Vararg{Int,N}) = ZeroPad{length(dim_in)}(Float64, dim_in, zp...)

ZeroPad(x::AbstractArray, zp::Vararg{Int}) = ZeroPad{ndims(x)}(eltype(x), size(x), zp...)

# Operators
function A_mul_B!{T}(y::AbstractVector{T}, L::ZeroPad{1}, b::AbstractVector{T}) 
	for i in eachindex(y)
		y[i] = i <= length(b) ? b[i] : 0.
	end
end

function Ac_mul_B!{T}(y::AbstractVector{T}, L::ZeroPad{1}, b::AbstractVector{T}) 
	for i in eachindex(y)
		y[i] = b[i] 
	end
end

function A_mul_B!{T}(y::AbstractArray{T,2}, L::ZeroPad{2}, b::AbstractArray{T,2}) 
	for l = 1:size(y,1), m = 1:size(y,2)
		y[l,m] = l <= size(b,1) && m <= size(b,2) ? b[l,m] : 0.
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,2}, L::ZeroPad{2}, b::AbstractArray{T,2}) 
	for l = 1:size(y,1), m = 1:size(y,2)
		y[l,m] = b[l,m]
	end
end

function A_mul_B!{T}(y::AbstractArray{T,3}, L::ZeroPad{3}, b::AbstractArray{T,3}) 
	for l = 1:size(y,1), m = 1:size(y,2), n = 1:size(y,3)
		y[l,m,n] = l <= size(b,1) && m <= size(b,2) && n <= size(b,3) ? b[l,m,n] : 0.
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,3}, L::ZeroPad{3}, b::AbstractArray{T,3}) 
	for l = 1:size(y,1), m = 1:size(y,2), n = 1:size(y,3)
		y[l,m,n] = b[l,m,n]
	end
end

fun_name(L::ZeroPad)       = "Zero Pad"
isGramDiagonal(L::ZeroPad)  = true
