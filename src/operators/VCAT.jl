import Base: vcat

immutable VCAT <: LinearOperator
	A::AbstractArray{LinearOperator}
	mid::AbstractArray
	function VCAT{T<:LinearOperator}(A::AbstractArray{T}, mid::AbstractArray)
		if any(size.(A[2:end], 2) .!= size(A[1], 2))
			throw(DimensionMismatch("operators must have the same codomain dimension!"))
		end
		if any(domainType.(A[2:end]) .!= domainType(A[1]))
			throw(DomainError())
		end
		new(A, mid)
	end
end
size(L::VCAT) = tuple(size.(L.A, 1)...), size(L.A[1],2)

# Constructors
VCAT(A::LinearOperator) = A

function VCAT(A::Vararg{LinearOperator})
	mid = zeros(domainType(A[1]),size(A[1],2))
	return VCAT([A...], mid)
end

## define `hcat` for convenience
vcat(L::Vararg{LinearOperator}) = VCAT(L...)
vcat(L::LinearOperator) = L

# Operators
function A_mul_B!{T1,T2<:AbstractArray}(y::AbstractArray{T2}, S::VCAT, b::AbstractArray{T1})
	for i = 1:length(S.A)
		A_mul_B!(y[i],S.A[i],b)
	end
end

function *(A::VCAT,b::AbstractArray)
	y = Array{AbstractArray,1}(length(A.A))
	for i = 1:length(A.A)
		C = codomainType(A.A[i])
		y[i] = Array{C}(size(A.A[i],1))
	end
	A_mul_B!(y,A,b)
	return y
end

# Transformations
transpose(L::VCAT) = HCAT((L.A.')[:],L.mid)

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

