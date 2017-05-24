export Compose

immutable Compose{N, M, L<:NTuple{N,Any}, T<:NTuple{M,Any}, C <: AbstractArray, D <: AbstractArray} <: LinearOperator
	A::L
	mid::T       # memory in the middle of the operators
end

# Constructors

function Compose(L1::LinearOperator, L2::LinearOperator)
	if size(L1,2) != size(L2,1)
		throw(DimensionMismatch("cannot compose operators"))
	end
	if domainType(L1) != codomainType(L2)
		throw(DomainError())
	end
	Compose( L1, L2, Array{domainType(L1)}(size(L2,1)) )
end

Compose(L1::LinearOperator,L2::LinearOperator,mid::AbstractArray) =
Compose( (L2,L1), (mid,))

Compose(L1::Compose,       L2::LinearOperator,mid::AbstractArray) =
Compose( (L2,L1.A...), (mid,L1.mid...))

Compose(L1::LinearOperator,L2::Compose,       mid::AbstractArray) =
Compose((L2.A...,L1), (L2.mid...,mid))

Compose(L1::Compose,       L2::Compose,       mid::AbstractArray) =
Compose((L2.A...,L1.A...), (L2.mid...,mid,L1.mid...))

Compose{N,M}(A::NTuple{N,Any},mid::NTuple{M,Any}) =
Compose{N,M,typeof(A),typeof(mid),
	Array{codomainType(A[end]),ndims(A[end],1)},
	Array{domainType(A[1]),ndims(A[1],2)}}(A,mid)

# Syntax (commented for now; does not belong here)

# *(L1::LinearOperator, L2::LinearOperator) = Compose(L1,L2)

# *{E<:Eye}(L1::E, L2::LinearOperator) = L2
# *{E<:Eye}(L1::LinearOperator, L2::E) = L1
# *{E1<:Eye,E2<:Eye}(L1::E1, L2::E2) = L1
#
# *{S<:Scale}(L1::S, L2::LinearOperator) = L1.coeff*(L1.A*L2)
# *{S<:Scale}(L1::LinearOperator, L2::S) = L2.coeff*(L1*L2.A)
# *{S1<:Scale,S2<:Scale}(L1::S1, L2::S2) = (L1.coeff*L2.coeff)*(L1.A*L2.A)
#
# # redefine .*
# Base.broadcast(::typeof(*), d::AbstractArray, L::LinearOperator) = DiagOp(codomainType(L), d)*L
# Base.broadcast(::typeof(*), d::AbstractArray, L::Scale)          = DiagOp(L.coeff*d)

# Mappings

@generated function A_mul_B!{N,M,T1,T2,C,D}(y::C, L::Compose{N,M,T1,T2,C,D},b::D)
	ex = :(A_mul_B!(L.mid[1],L.A[1],b))
	for i = 2:M
		ex = quote
			$ex
			A_mul_B!(L.mid[$i],L.A[$i], L.mid[$i-1])
		end
	end
	ex = quote
		$ex
		A_mul_B!(y,L.A[N], L.mid[M])
		return y
	end
end

@generated function Ac_mul_B!{N,M,T1,T2,C,D}(y::D, L::Compose{N,M,T1,T2,C,D},b::C)
	ex = :(Ac_mul_B!(L.mid[M],L.A[N],b))
	for i = M:-1:2
		ex = quote
			$ex
			Ac_mul_B!(L.mid[$i-1],L.A[$i], L.mid[$i])
		end
	end
	ex = quote
		$ex
		Ac_mul_B!(y,L.A[1], L.mid[1])
		return y
	end
end

# Properties

size(L::Compose) = ( size(L.A[end],1), size(L.A[1],2) )

fun_name(L::Compose) = length(L.A) == 2 ? fun_name(L.A[2])*" * "*fun_name(L.A[1]) : "Composition"

domainType(L::Compose)   = domainType(L.A[1])
codomainType(L::Compose) = codomainType(L.A[end])

is_diagonal(L::Compose) = all(is_diagonal.(L.A))
is_invertible(L::Compose) = all(is_invertible.(L.A))
