import Base: vcat

immutable VCAT{N, C<:NTuple{N,Any}, D<:AbstractArray, L<:NTuple{N,Any}} <: LinearOperator
	A::L
	mid::D
end
size(L::VCAT) = size.(L.A, 1), size(L.A[1],2)

# Constructors
function VCAT{N,D<:AbstractArray}(A::NTuple{N,Any}, mid::D)
	if any([size(A[1],2) != size(a,2) for a in A]) 
		throw(DimensionMismatch("operators must have the same codomain dimension!"))
	end
	if any(domainType.(A[2:end]) .!= domainType(A[1]))
		throw(DomainError())
	end
	codomType = codomainType.(A) 
	C = Tuple{[Array{codomType[i],ndims(A[i],1)} for i in eachindex(codomType)]...}
	VCAT{N,C,D,typeof(A)}(A, mid)
end

VCAT(A::LinearOperator) = A

function VCAT(A::Vararg{LinearOperator})
	mid = zeros(domainType(A[1]),size(A[1],2))
	return VCAT(A, mid)
end

## define `hcat` for convenience
vcat(L::Vararg{LinearOperator}) = VCAT(L...)
vcat(L::LinearOperator) = L

# Operators
function A_mul_B!{N,C,D,L}(y::C, V::VCAT{N,C,D,L}, b::D)
	for i = 1:length(V.A)
		A_mul_B!(y[i],V.A[i],b)
	end
end

function *(V::VCAT,b::AbstractArray)
	y = zeros.(codomainType.(V.A),size.(V.A,1))
	A_mul_B!(y,V,b)
	return y
end

# Transformations
function transpose{N,C,D,L}(V::VCAT{N,C,D,L}) 
	At = transpose.(V.A)
	HCAT{N,D,C,typeof(At)}(At,V.mid)
end

# Properties
fun_name(L::VCAT) = length(L.A) == 2 ? "["fun_name(L.A[1])*"; "*fun_name(L.A[2])*"]"  :"VCAT operator"

# Utils

(+)(L1::VCAT, L2::VCAT) = VCAT(L1.A.+ L2.A, L1.mid)
(-)(L1::VCAT, L2::VCAT) = VCAT(L1.A.-L2.A, L1.mid)

function fun_codomain(L::VCAT)
	str = ""
	for i in eachindex(L.A) 
		str *= fun_codomain(L.A[i])
		i != length(L.A) && (str *= ", ")
	end
	return str
end

  domainType(L::VCAT) = domainType(L.A[1])
codomainType(L::VCAT) = codomainType.(L.A)

