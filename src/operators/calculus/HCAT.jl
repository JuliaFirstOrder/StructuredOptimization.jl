export HCAT

immutable HCAT{M, N,
	       C <: Union{NTuple{M,AbstractArray}, AbstractArray},
	       D <: NTuple{N,AbstractArray},
	       L <: NTuple{N,LinearOperator}} <: LinearOperator
	A::L
	mid::C
end

# Constructors

function HCAT{N, C<:Union{Tuple,AbstractArray}, L <:NTuple{N,LinearOperator}}(A::L, mid::C, M::Int)
	if any([size(A[1],1) != size(a,1) for a in A])
		throw(DimensionMismatch("operators must have the same codomain dimension!"))
	end
	if any([codomainType(A[1]) != codomainType(a) for a in A])
		throw(DomainError())
	end
	domType = domainType.(A)
	D = Tuple{[Array{domType[i],ndims(A[i],2)} for i in eachindex(domType)]...}
	HCAT{M,N,C,D,L}(A, mid)
end

HCAT(A::LinearOperator) = A

function HCAT(A::Vararg{LinearOperator})
	s = size(A[1],1)
	t = codomainType(A[1])
	mid,M  = create_mid(t,s)
	return HCAT(A, mid, M)
end

create_mid{N}(t::NTuple{N,DataType},s::NTuple{N,NTuple}) = zeros.(t,s), N
create_mid{N}(t::Type,s::NTuple{N,Int}) = zeros(t,s), 1

# Syntax (commented for now; does not belong here)

# import Base: hcat
# hcat(L::Vararg{LinearOperator}) = HCAT(L...)
# hcat(L::LinearOperator) = L
# (+)(L1::HCAT, L2::HCAT) = HCAT(L1.A.+ L2.A, L1.mid)
# (-)(L1::HCAT, L2::HCAT) = HCAT(L1.A.-L2.A, L1.mid)

# Mappings

@generated function A_mul_B!{M,N,C,D,L}(y::C, S::HCAT{M,N,C,D,L}, b::D)
	ex = :(A_mul_B!(y, S.A[1], b[1]))
	for i = 2:N
		ex = quote
			$ex
			A_mul_B!(S.mid, S.A[$i], b[$i])
		end

		if C <: AbstractArray
			ex = :($ex; y .+= S.mid)
		else
			for ii = 1:M
				ex = :($ex; y[$ii] .+= S.mid[$ii])
			end
		end
	end
	ex = quote
		$ex
		return y
	end
end

@generated function Ac_mul_B!{M,N,C,D,L}(y::D, H::HCAT{M,N,C,D,L}, b::C)
	ex = :()
	for i = 1:N
		ex = :($ex; Ac_mul_B!(y[$i],H.A[$i],b))
	end
	ex = quote
		$ex
		return y
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

fun_codomain(L::HCAT) = fun_codomain(L.A[1])

domainType(L::HCAT) = domainType.(L.A)
codomainType(L::HCAT) = codomainType.(L.A[1])

is_gram_diagonal(L::HCAT) = all(is_gram_diagonal.(L.A))
is_full_row_rank(L::HCAT) = any(is_full_row_rank.(L.A))
