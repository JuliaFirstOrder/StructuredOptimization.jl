export Reshape

immutable Reshape{N,M,L<:LinearOperator} <: LinearOperator
	dim_out::NTuple{N,Int}
	dim_in::NTuple{M,Int}
	A::L
end

size(L::Reshape) = (L.dim_out, L.dim_in)

# Constructors
Reshape(L::LinearOperator, dim_out...) = 
Reshape{length(dim_out),length(size(L,2)),typeof(L)}( dim_out, size(L,2), L)
Reshape(L::Scale, dim_out...) = 
L.coeff*(Reshape( L.A, dim_out...))

# Operators
function A_mul_B!(y::AbstractArray,L::Reshape,b::AbstractArray)
	y_res = reshape(y,size(L.A,1))
	b_res = reshape(b,size(L.A,2))
	A_mul_B!(y_res, L.A, b_res)
end

function Ac_mul_B!(y::AbstractArray,L::Reshape,b::AbstractArray)
	y_res = reshape(y,size(L.A,2))
	b_res = reshape(b,size(L.A,1))
	Ac_mul_B!(y_res, L.A, b_res)
end

# Transformations

# Properties
  domainType(  L::Reshape) =   domainType(L.A)
codomainType(  L::Reshape) = codomainType(L.A)
isEye(         L::Reshape) = isEye(L.A) 
isDiagonal(    L::Reshape) = isDiagonal(L.A) 
isGramDiagonal(L::Reshape) = isGramDiagonal(L.A)
isInvertible(  L::Reshape) = isInvertible(L.A)

fun_name(L::Reshape) = "Reshaped "*fun_name(L.A)

