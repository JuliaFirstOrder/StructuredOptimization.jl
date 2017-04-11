export Reshape

immutable Reshape <: LinearOperator
	dim_out::Tuple
	dim_in::Tuple
	A::LinearOperator
end

size(L::Reshape) = (L.dim_out, L.dim_in)

# Constructors
Reshape(L::LinearOperator, dim_out...) = Reshape( dim_out, size(L,2), L)
Reshape(L::Scale, dim_out...) = L.coeff*(Reshape( dim_out, size(L,2), L.A))

# Operators
function A_mul_B!(y::AbstractArray,L::Reshape,b::AbstractArray)
	y_res = reshape(y,size(L.A,1))
	b_res = reshape(b,size(L.A,2))
	A_mul_B!(y_res, L.A, b_res)
end

# Transformations
transpose(L::Reshape) = Reshape(L.dim_in,L.dim_out, L.A')

# Properties
  domainType(  L::Reshape) =   domainType(L.A)
codomainType(  L::Reshape) = codomainType(L.A)
isEye(         L::Reshape) = isEye(L.A) 
isDiagonal(    L::Reshape) = isDiagonal(L.A) 
isGramDiagonal(L::Reshape) = isGramDiagonal(L.A)
isInvertible(  L::Reshape) = isInvertible(L.A)

fun_name(L::Reshape) = "Reshaped "*fun_name(L.A)

