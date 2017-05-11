import Base: +, -

immutable Sum{N,C<:AbstractArray,D<:AbstractArray,T<:NTuple{N,Any}} <: LinearOperator
	A::T
	mid::C
end
	
function Sum{N}(A::NTuple{N,Any})
	if any(  [size(a) != size(A[1]) for a in A]   )
		throw(DimensionMismatch("cannot sum operator of different sizes"))
	end
	if any(codomainType.(A) .!= codomainType(A[1])) || 
		any(  domainType.(A) .!=   domainType(A[1]))
		throw(DomainError())
	end
	mid = Array{codomainType(A[1])}(size(A[1],1))
	Sum{N,typeof(mid),Array{domainType(A[1]),length(size(A[1],2))},typeof(A)}(A, mid)
end

size(L::Sum) = size(L.A[1])

# Constructors
-(L::Sum) = Sum((-).(L.A), L.mid)

+(L1::LinearOperator, L2::LinearOperator) = Sum((L1,  L2 ))
-(L1::LinearOperator, L2::LinearOperator) = Sum((L1, -L2 ))

+(L1::LinearOperator, L2::Sum)     = Sum((L1,L2.A...))
-(L1::LinearOperator, L2::Sum)     = Sum((L1,((-).(L2.A))...))

+(L1::Sum, L2::LinearOperator) = L2+L1
-(L1::Sum, L2::LinearOperator) = Sum((L1.A..., -L2))

# Operators
@generated function A_mul_B!{N,C,D}(y::C, S::Sum{N,C,D}, b::D)
	ex = :(A_mul_B!(y,S.A[1],b))
	for i = 2:N
		ex = quote 
			$ex
			A_mul_B!(S.mid,S.A[$i],b)
			y .+= S.mid
		end
	end
	ex = quote
		$ex
		return y
	end
end

# Transformations
transpose{N,C,D,T}(S::Sum{N,C,D,T}) = 
Sum{N,D,C,typeof(transpose.(S.A))}(transpose.(S.A),Array{domainType(S.A[1])}(size(S,2)))

# Properties

domainType(L::Sum) =     domainType(L.A[1])
codomainType(L::Sum) = codomainType(L.A[1])

fun_name(S::Sum) = 
length(S.A) == 2 ? fun_name(S.A[1])" + "fun_name(S.A[2]) : "Sum of linear operators"
