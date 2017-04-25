import Base: hcat

immutable HCAT <: LinearOperator
	A::AbstractArray{LinearOperator}
	mid::AbstractArray
	function HCAT{T<:LinearOperator}(A::AbstractArray{T}, mid::AbstractArray)
		if any(size.(A[2:end], 1) .!= size(A[1], 1))
			throw(DimensionMismatch("operators must have the same codomain dimension!"))
		end
		if any(codomainType.(A[2:end]) .!= codomainType(A[1]))
			throw(DomainError())
		end
		new(A, mid)
	end
end
size(L::HCAT) = size(L.A[1],1), tuple(size.(L.A, 2)...)

# Constructors
HCAT(A::LinearOperator) = A

function HCAT(A::Vararg{LinearOperator})
	mid = zeros(codomainType(A[1]),size(A[1],1))
	return HCAT([A...], mid)
end

## define `hcat` for convenience
hcat(L::Vararg{LinearOperator}) = HCAT(L...)
hcat(L::LinearOperator) = L

# Operators
function A_mul_B!{T1,T2<:AbstractArray}(y::AbstractArray{T1}, S::HCAT, b::AbstractArray{T2})
	A_mul_B!(y, S.A[1], b[1])
	for i = 2:length(S.A)
		A_mul_B!(S.mid, S.A[i], b[i])
		y .= (+).(y, S.mid)
	end
end

# Transformations
transpose(L::HCAT) = VCAT((L.A.')[:],L.mid)

# Properties
fun_name(L::HCAT) = length(L.A) == 2 ? "["fun_name(L.A[1])*", "*fun_name(L.A[2])*"]"  : "HCAT operator"

# Utils

(+)(L1::HCAT, L2::HCAT) = HCAT(L1.A.+ L2.A, L1.mid)
(-)(L1::HCAT, L2::HCAT) = HCAT(L1.A.-L2.A, L1.mid)

function fun_domain(L::HCAT)
	str = ""
	for i in eachindex(L.A) 
		str *= fun_domain(L.A[i])
		i != length(L.A) && (str *= ", ")
	end
	return str
end

  domainType(L::HCAT) = domainType.(L.A)
codomainType(L::HCAT) = codomainType(L.A[1])


# import Base: copy, sort
#
# copy(A::HCAT) = HCAT(copy(A.A), A.mid)
#
# function sort(A::HCAT,p::Array)
# 	H = A.A[p]
# 	return HCAT(H,A.mid)
# end
#
