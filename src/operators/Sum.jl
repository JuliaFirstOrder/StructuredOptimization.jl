import Base: +, -

immutable Sum <: LinearOperator
	A::Vector{LinearOperator}
	mid::AbstractArray
	function Sum(A, mid)
		if any(size.(A[2:end]) .!= size(A[1]))
			throw(DimensionMismatch("cannot sum operator of different sizes"))
		end
		if any(codomainType.(A[2:end]) .!= codomainType(A[1])) || 
		   any(  domainType.(A[2:end]) .!=   domainType(A[1]))
			throw(DomainError())
		end
		new(A, mid)
	end
end

size(L::Sum) = size(L.A[1])

# Constructors
-(L::Sum) = Sum((-).(L.A), L.mid)

+(L1::LinearOperator, L2::LinearOperator) = Sum([L1,  L2 ], Array{codomainType(L1)}(size(L1,1)))
-(L1::LinearOperator, L2::LinearOperator) = Sum([L1,(-L2)], Array{codomainType(L1)}(size(L1,1)))

+(L1::LinearOperator, L2::Sum)     = Sum([L1,L2.A...],L2.mid)
-(L1::LinearOperator, L2::Sum)     = Sum([L1,(-L2.A)...],L2.mid)

+(L1::Sum, L2::LinearOperator) = L2+L1
-(L1::Sum, L2::LinearOperator) = Sum([L1.A...,(-L2)],L1.mid)

# Operators
function A_mul_B!(y::AbstractArray,S::Sum,b::AbstractArray)
	A_mul_B!(y,S.A[1],b)
	for i = 2:length(S.A)
		A_mul_B!(S.mid,S.A[i],b)
		y .= (+).(y,S.mid)
	end
end

# Transformations
transpose(S::Sum) = Sum((S.A.')[:],Array{domainType(S.A[1])}(size(S,2)))

# Properties

domainType(L::Sum) =     domainType(L.A[1])
codomainType(L::Sum) = codomainType(L.A[1])

fun_name(S::Sum) = 
length(S.A) == 2 ? fun_name(S.A[1])" + "fun_name(S.A[2]) : "Sum of linear operators"
