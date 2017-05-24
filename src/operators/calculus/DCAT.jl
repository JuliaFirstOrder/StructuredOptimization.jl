export DCAT

immutable DCAT{N,
	       C <: NTuple{N,AbstractArray},
	       D <: NTuple{N,AbstractArray},
	       L<:NTuple{N,LinearOperator}} <: LinearOperator
	A::L
end

# Constructors
DCAT(A::LinearOperator) = A

function DCAT(A::Vararg{LinearOperator})
	domType, codomType = domainType.(A), codomainType.(A)
	C = Tuple{[Array{codomType[i],ndims(A[i],1)} for i in eachindex(codomType)]...}
	D = Tuple{[Array{  domType[i],ndims(A[i],2)} for i in eachindex(  domType)]...}
	return DCAT{length(A), C, D, typeof(A)}(A)
end

# Mappings
A_mul_B!{N,C,D,L}(y::C, S::DCAT{N,C,D,L}, b::D) = A_mul_B!.(y,S.A,b)

Ac_mul_B!{N,C,D,L}(y::D, S::DCAT{N,C,D,L}, b::C) = Ac_mul_B!.(y,S.A,b)

# Properties
size(L::DCAT) = size.(L.A,1), size.(L.A, 2)

fun_name(L::DCAT) = length(L.A) == 2 ? "["fun_name(L.A[1])*", 0; 0 "*fun_name(L.A[2])*"]" :
"DCAT operator"

function fun_domain(L::DCAT)
	str = ""
	for i in eachindex(L.A)
		str *= fun_domain(L.A[i])
		i != length(L.A) && (str *= ", ")
	end
	return str
end

function fun_codomain(L::DCAT)
	str = ""
	for i in eachindex(L.A)
		str *= fun_codomain(L.A[i])
		i != length(L.A) && (str *= ", ")
	end
	return str
end

domainType(L::DCAT)   = domainType.(L.A)
codomainType(L::DCAT) = codomainType.(L.A)

is_diagonal(L::DCAT) = all(is_diagonal.(L.A))
is_gram_diagonal(L::DCAT) = all(is_gram_diagonal.(L.A))
is_invertible(L::DCAT) = all(is_invertible.(L.A))
is_full_row_rank(L::DCAT) = all(is_full_row_rank.(L.A))
is_full_column_rank(L::DCAT) = all(is_full_column_rank.(L.A))
