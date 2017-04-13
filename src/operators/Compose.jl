immutable Compose <: LinearOperator
	A::Vector{LinearOperator}
	mid::Vector{AbstractArray}       # memory in the middle of the operators

	function Compose(L1::LinearOperator, L2::LinearOperator)
		if size(L1,2) != size(L2,1)
			throw(DimensionMismatch("cannot compose operators"))
		end
		if domainType(L1) != codomainType(L2)
			throw(DomainError())
		end
		Compose(L1,L2,Array{domainType(L1)}(size(L2,1)))
	end

	Compose(L1::LinearOperator,L2::LinearOperator,mid::AbstractArray) = new([L2,L1], [mid])

	Compose(L1::Compose,       L2::LinearOperator,mid::AbstractArray) = 
	new([L2,L1.A...], [mid,L1.mid...])

	Compose(L1::LinearOperator,L2::Compose,       mid::AbstractArray) = 
	new([L2.A...,L1], [L2.mid...,mid])

	Compose(L1::Compose,       L2::Compose,       mid::AbstractArray) = 
	new([L2.A...,L1.A...], [L2.mid...,mid,L1.mid...])

	Compose(A,mid) = new(A,mid)
end

size(L::Compose) = ( size(L.A[end],1), size(L.A[1],2) )

# Constructors

*(L1::LinearOperator, L2::LinearOperator) = Compose(L1,L2)

# Special cases

*{E<:IdentityOperator}(L1::E, L2::LinearOperator) = L2
*{E<:IdentityOperator}(L1::LinearOperator, L2::E) = L1
*{E1<:IdentityOperator,E2<:IdentityOperator}(L1::E1, L2::E2) = L1

*{S<:Scale}(L1::S, L2::LinearOperator) = L1.coeff*(L1.A*L2)
*{S<:Scale}(L1::LinearOperator, L2::S) = L2.coeff*(L1*L2.A)
*{S1<:Scale,S2<:Scale}(L1::S1, L2::S2) = (L1.coeff*L2.coeff)*(L1.A*L2.A)

.*(d::AbstractArray, L2::LinearOperator) = DiagOp(codomainType(L2), d)*L2
.*(d::AbstractArray, L2::Scale         ) = DiagOp(L2.coeff*d)

# Operators
function A_mul_B!(y::AbstractArray,L::Compose,b::AbstractArray)
	A_mul_B!(L.mid[1],L.A[1],b)
	for i = 2:length(L.A)-1
		A_mul_B!(L.mid[i],L.A[i], L.mid[i-1])
	end
	A_mul_B!(y,L.A[length(L.A)], L.mid[length(L.A)-1])
end

# Transformations
function transpose(L::Compose)
	Compose(flipdim((L.A.')[:],1),flipdim(L.mid,1))
end

# Properties
fun_name(L::Compose) = length(L.A) == 2 ? fun_name(L.A[2])*" * "*fun_name(L.A[1]) : "Nested Linear Operator"

domainType(L::Compose)   = domainType(L.A[1])
codomainType(L::Compose) = codomainType(L.A[end])

