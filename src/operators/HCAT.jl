import Base: hcat

immutable HCAT{N, C<:AbstractArray, D<:NTuple{N,Any}, L<:NTuple{N,Any}} <: LinearOperator
	A::L
	mid::C
end
size(L::HCAT) = size(L.A[1],1), size.(L.A, 2)

# Constructors
function HCAT{N,C<:AbstractArray}(A::NTuple{N,Any}, mid::C)
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
	mid = zeros(codomainType(A[1]),size(A[1],1))
	return HCAT(A, mid)
end

## define `hcat` for convenience
hcat(L::Vararg{LinearOperator}) = HCAT(L...)
hcat(L::LinearOperator) = L

# Operators
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

#function A_mul_B!{N,C,D,L}(y::C, S::HCAT{N,C,D,L}, b::D)
#	A_mul_B!(y, S.A[1], b[1])
#	for i = 2:length(S.A)
#		A_mul_B!(S.mid, S.A[i], b[i])
#		y .= (+).(y, S.mid)
#	end
#end

function (*){N,C,D,L}(H::HCAT{N,C,D,L},b::D)
	y = zeros(codomainType(H),size(H,1))
	A_mul_B!(y,H,b)
	return y
end

# Transformations
function transpose{N,C,D,L}(H::HCAT{N,C,D,L}) 
	At = transpose.(H.A)
	VCAT{N,D,C,typeof(At)}(At,H.mid)
end

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

