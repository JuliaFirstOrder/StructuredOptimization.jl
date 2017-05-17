export Sum

immutable Sum{N, C<:AbstractArray, D<:AbstractArray, T<:NTuple{N,Any}} <: LinearOperator
	A::T
	midC::C
	midD::D
end

# Constructors

function Sum{N}(A::NTuple{N,Any}, midC, midD)
	if any([size(a) != size(A[1]) for a in A])
		throw(DimensionMismatch("cannot sum operator of different sizes"))
	end
	if any(codomainType.(A) .!= codomainType(A[1])) ||
		any(domainType.(A) .!= domainType(A[1]))
		throw(DomainError())
	end
	Sum{N, typeof(midC), typeof(midD), typeof(A)}(A, midC, midD)
end

Sum(A::LinearOperator) = A

function Sum(A::Vararg{LinearOperator})
	midC = Array{codomainType(A[1])}(size(A[1], 1))
	midD = Array{domainType(A[1])}(size(A[1], 2))
	return Sum(A, midC, midD)
end

# Syntax (commented for now; does not belong here)

# import Base: +, -
# -(L::Sum) = Sum((-).(L.A), L.midC, L.midD)
# +(L1::LinearOperator, L2::LinearOperator) = Sum((L1,  L2 ))
# -(L1::LinearOperator, L2::LinearOperator) = Sum((L1, -L2 ))
# +(L1::LinearOperator, L2::Sum)     = Sum((L1,L2.A...))
# -(L1::LinearOperator, L2::Sum)     = Sum((L1,((-).(L2.A))...))
# +(L1::Sum, L2::LinearOperator) = L2+L1
# -(L1::Sum, L2::LinearOperator) = Sum((L1.A..., -L2))

# Mappings

@generated function A_mul_B!{N,C,D}(y::C, S::Sum{N,C,D}, b::D)
	ex = :(A_mul_B!(y,S.A[1],b))
	for i = 2:N
		ex = quote
			$ex
			A_mul_B!(S.midC,S.A[$i],b)
			y .+= S.midC
		end
	end
	ex = quote
		$ex
		return y
	end
end

@generated function Ac_mul_B!{N,C,D}(y::D, S::Sum{N,C,D}, b::C)
	ex = :(Ac_mul_B!(y,S.A[1],b))
	for i = 2:N
		ex = quote
			$ex
			Ac_mul_B!(S.midD,S.A[$i],b)
			y .+= S.midD
		end
	end
	ex = quote
		$ex
		return y
	end
end

# Properties

size(L::Sum) = size(L.A[1])

domainType(L::Sum) = domainType(L.A[1])
codomainType(L::Sum) = codomainType(L.A[1])

fun_name(S::Sum) =
length(S.A) == 2 ? fun_name(S.A[1])" + "fun_name(S.A[2]) : "Sum of linear operators"
