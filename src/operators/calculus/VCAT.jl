export VCAT

immutable VCAT{N, C <: NTuple{N,Any}, D <: AbstractArray, L <: NTuple{N,Any}} <: LinearOperator
	A::L
	mid::D
end

# Constructors

function VCAT{N, D <: AbstractArray}(A::NTuple{N,Any}, mid::D)
	if any([size(A[1],2) != size(a,2) for a in A])
		throw(DimensionMismatch("operators must have the same codomain dimension!"))
	end
	if any(domainType.(A[2:end]) .!= domainType(A[1]))
		throw(DomainError())
	end
	codomType = codomainType.(A)
	C = Tuple{[Array{codomType[i],ndims(A[i],1)} for i in eachindex(codomType)]...}
	VCAT{N, C, D, typeof(A)}(A, mid)
end

VCAT(A::LinearOperator) = A

function VCAT(A::Vararg{LinearOperator})
	mid = zeros(domainType(A[1]), size(A[1], 2))
	return VCAT(A, mid)
end

# Syntax (commented for now; does not belong here)

# import Base: vcat
# vcat(L::Vararg{LinearOperator}) = VCAT(L...)
# vcat(L::LinearOperator) = L
# (+)(L1::VCAT, L2::VCAT) = VCAT(L1.A.+ L2.A, L1.mid)
# (-)(L1::VCAT, L2::VCAT) = VCAT(L1.A.-L2.A, L1.mid)

# Mappings

function A_mul_B!{N,C,D,L}(y::C, V::VCAT{N,C,D,L}, b::D)
	for i = 1:length(V.A)
		A_mul_B!(y[i],V.A[i],b)
	end
end

@generated function Ac_mul_B!{N,C,D,L}(y::D, S::VCAT{N,C,D,L}, b::C)
	ex = :(Ac_mul_B!(y, S.A[1], b[1]))
	for i = 2:N
		ex = quote
			$ex
			Ac_mul_B!(S.mid, S.A[$i], b[$i])
			y .+= S.mid
		end
	end
	ex = quote
		$ex
		return y
	end

end

# Properties

size(L::VCAT) = size.(L.A, 1), size(L.A[1],2)

fun_name(L::VCAT) = length(L.A) == 2 ? "["fun_name(L.A[1])*"; "*fun_name(L.A[2])*"]"  :"VCAT operator"

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
