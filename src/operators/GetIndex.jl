import Base: getindex
export GetIndex 

immutable GetIndex{N,M,T<:Tuple} <: LinearOperator
	domainType::Type
	dim_out::NTuple{N,Int}
	dim_in::NTuple{M,Int}
	idx::T
end

size(L::GetIndex) = (L.dim_out,L.dim_in)

# Constructors

function GetIndex{M,T<:Tuple}(domainType::Type,dim_in::NTuple{M,Int},idx::T)
	dim_out = get_dim_out(dim_in,idx...)
	GetIndex{length(dim_out),M,T}(domainType,dim_out,dim_in,idx)
end

GetIndex(dim_in::Tuple, idx::Tuple) = GetIndex(Float64, dim_in, idx)
GetIndex(x::AbstractArray, idx::Tuple) = GetIndex(eltype(x), size(x), idx)

getindex(A::LinearOperator,idx...) = GetIndex(codomainType(A),size(A,1),idx)*A 

# Operators
function A_mul_B!{T1,N,M,T2}(y::Array{T1,N},L::GetIndex{N,M,T2},b::Array{T1,M})
	y .= view(b,L.idx...)
end

function Ac_mul_B!{T1,N,M,T2}(y::Array{T1,M},L::GetIndex{N,M,T2},b::AbstractArray{T1,N})
	y .= 0.
	setindex!(y,b,L.idx...)
end

# Properties
fun_name(L::GetIndex) = "Get Index"
isGramDiagonal(L::GetIndex) = true

# Utils
get_idx(L::GetIndex) = L.idx

function get_dim_out(dim,args...)
	if length(args) != 1
		dim2 = [dim...]
		for i = 1:length(args)
			if args[i] != Colon() dim2[i] = length(args[i]) end
		end
		return tuple(dim2...)
	else
		return tuple(length(args[1]))
	end
end
