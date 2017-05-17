export HCAT

immutable HCAT{N, C <: AbstractArray, D <: NTuple{N,Any}, L <: NTuple{N,Any}} <: LinearOperator
	A::L
	mid::C
end

# Constructors

function HCAT{N, C <: AbstractArray}(A::NTuple{N,Any}, mid::C)
	if any([size(A[1],1) != size(a,1) for a in A])
		throw(DimensionMismatch("operators must have the same codomain dimension!"))
	end
	if any(codomainType.(A[2:end]) .!= codomainType(A[1]))
		throw(DomainError())
	end
	domType = domainType.(A)
	D = Tuple{[Array{domType[i],ndims(A[i],2)} for i in eachindex(domType)]...}
	HCAT{N,C,D,typeof(A)}(A, mid)
end

HCAT(A::LinearOperator) = A

function HCAT(A::Vararg{LinearOperator})
	mid = zeros(codomainType(A[1]), size(A[1], 1))
	return HCAT(A, mid)
end

# Syntax (commented for now; does not belong here)

# import Base: hcat
# hcat(L::Vararg{LinearOperator}) = HCAT(L...)
# hcat(L::LinearOperator) = L
# (+)(L1::HCAT, L2::HCAT) = HCAT(L1.A.+ L2.A, L1.mid)
# (-)(L1::HCAT, L2::HCAT) = HCAT(L1.A.-L2.A, L1.mid)

# Mappings

@generated function A_mul_B!{N,C,D,L}(y::C, S::HCAT{N,C,D,L}, b::D)
	ex = :(A_mul_B!(y, S.A[1], b[1]))
	for i = 2:N
		ex = quote
			$ex
			A_mul_B!(S.mid, S.A[$i], b[$i])
			y .+= S.mid
		end
	end
	ex = quote
		$ex
		return y
	end
end

function Ac_mul_B!{N,C,D,L}(y::D, H::HCAT{N,C,D,L}, b::C)
	for i = 1:length(H.A)
		Ac_mul_B!(y[i],H.A[i],b)
	end
end

# Properties

size(L::HCAT) = size(L.A[1],1), size.(L.A, 2)

fun_name(L::HCAT) = length(L.A) == 2 ? "["fun_name(L.A[1])*", "*fun_name(L.A[2])*"]" : "HCAT operator"

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
