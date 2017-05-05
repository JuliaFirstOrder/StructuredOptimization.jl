import Base: getindex
export GetIndex 

immutable GetIndex{No,Ni} <: LinearOperator
	domainType::Type
	dim_in::NTuple{Ni,Int}
	dim_out::NTuple{No,Int}
	idx::Tuple
end

function GetIndex{N}(domainType::Type,dim_in::NTuple{N,Int},idx)
	dim_out = get_dim_out(dim_in,idx...)
	GetIndex{length(dim_out),N}(domainType,dim_in,dim_out,idx)
end

size(L::GetIndex) = (L.dim_out,L.dim_in)

# Constructors
GetIndex(dim_in::Tuple, idx::Tuple) = GetIndex(Float64, dim_in, idx)
GetIndex(x::AbstractArray, idx::Tuple) = GetIndex(eltype(x), size(x), idx)

getindex(A::LinearOperator,idx...) = GetIndex(codomainType(A),size(A,1),idx)*A 

# Operators
function A_mul_B!{T,No,Ni}(y::AbstractArray{T,No},L::GetIndex{No,Ni},b::AbstractArray{T,Ni})
	y .= view(b,L.idx...)
end

function Ac_mul_B!{T,No,Ni}(y::AbstractArray{T,Ni},L::GetIndex{No,Ni},b::AbstractArray{T,No})
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
