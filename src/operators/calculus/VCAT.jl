export VCAT

immutable VCAT{M, N, 
	       C <: NTuple{M,AbstractArray}, 
	       D <: Union{NTuple{N,AbstractArray}, AbstractArray}, 
	       L <: NTuple{M,LinearOperator}} <: LinearOperator
	A::L
	mid::D
end

# Constructors

function VCAT{M, D<:Union{Tuple,AbstractArray}, L<:NTuple{M,LinearOperator}}(A::L, mid::D, N::Int)
	if any([size(A[1],2) != size(a,2) for a in A])
		throw(DimensionMismatch("operators must have the same codomain dimension!"))
	end
	if any([domainType(A[1]) != domainType(a) for a in A])
		throw(DomainError())
	end
	codomType = codomainType.(A)
	C = Tuple{[Array{codomType[i],ndims(A[i],1)} for i in eachindex(codomType)]...}
	VCAT{M,N,C,D,L}(A, mid)
end

VCAT(A::LinearOperator) = A

function VCAT(A::Vararg{LinearOperator})
	s = size(A[1],2)
	t = domainType(A[1])
	mid,N  = create_mid(t,s)
	return VCAT(A, mid, N)
end

# Syntax (commented for now; does not belong here)

# import Base: vcat
# vcat(L::Vararg{LinearOperator}) = VCAT(L...)
# vcat(L::LinearOperator) = L
# (+)(L1::VCAT, L2::VCAT) = VCAT(L1.A.+ L2.A, L1.mid)
# (-)(L1::VCAT, L2::VCAT) = VCAT(L1.A.-L2.A, L1.mid)

# Mappings

@generated function A_mul_B!{M,N,C,D,L}(y::C, H::VCAT{M,N,C,D,L}, b::D)
	ex = :()
	for i = 1:M
		ex = :($ex; A_mul_B!(y[$i],H.A[$i],b))
	end
	ex = quote
		$ex
		return y
	end
end

@generated function Ac_mul_B!{M,N,C,D,L}(y::D, S::VCAT{M,N,C,D,L}, b::C)
	ex = :(Ac_mul_B!(y, S.A[1], b[1]))
	for i = 2:M
		ex = quote
			$ex
			Ac_mul_B!(S.mid, S.A[$i], b[$i])
		end
	
		if D <: AbstractArray
			ex = :($ex; y .+= S.mid)
		else
			for ii = 1:N
				ex = :($ex; y[$ii] .+= S.mid[$ii])
			end
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

fun_domain(L::VCAT) = fun_domain(L.A[1])

  domainType(L::VCAT) = domainType(L.A[1])
codomainType(L::VCAT) = codomainType.(L.A)
