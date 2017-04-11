export MatrixOp

immutable MatrixOp <: LinearOperator
	domainType::Type
	dim_in::Tuple
	A::AbstractMatrix

	MatrixOp(domainType,dim_in,A) = 
	new(domainType, dim_in, domainType <: Complex ? complex(A) : A)
end

size(L::MatrixOp) = length(L.dim_in) == 1 ? ((size(L.A,1),), L.dim_in) : 
((size(L.A,1), L.dim_in[2]), L.dim_in)

#Constructors

MatrixOp(A::AbstractMatrix) = MatrixOp(eltype(A), (size(A, 2),), A)

MatrixOp(A::AbstractMatrix, dim_in::Tuple) = MatrixOp(eltype(A), dim_in, A)

MatrixOp(A::AbstractMatrix, dim_in::Vararg{Int} ) = MatrixOp(A,dim_in) 

MatrixOp(domainType::Type, A::AbstractMatrix, dim_in::Vararg{Int} ) = MatrixOp(domainType,dim_in,A) 

MatrixOp(x::AbstractArray, A::AbstractMatrix) = MatrixOp(eltype(x),size(x),A)

# Operators
A_mul_B!{T}(y::AbstractArray{T},  A::MatrixOp, b::AbstractArray{T}) = A_mul_B!(y, A.A, b)
Ac_mul_B!{T}(y::AbstractArray{T}, A::MatrixOp, b::AbstractArray{T}) = Ac_mul_B!(y, A.A, b)

# Properties
fun_name(A::MatrixOp)  = "Matrix operator"

