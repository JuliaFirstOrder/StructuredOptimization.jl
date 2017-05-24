export MatrixOp

immutable MatrixOp{M <: AbstractMatrix, T} <: LinearOperator
	# domainType::Type
	# dim_in::NTuple{N,Int}
	A::M
	n_columns_input::Integer
end

# Constructors

# MatrixOp(domainType,dim_in,A) =
# MatrixOp(domainType, dim_in, domainType <: Complex ? complex(A) : A)

MatrixOp{M <: AbstractMatrix}(A::M) = MatrixOp{M, eltype(A)}(A, 1)
MatrixOp{M <: AbstractMatrix}(A::M, T::Type) = MatrixOp{M, T}(A, 1)
MatrixOp{M <: AbstractMatrix}(A::M, n::Integer) = MatrixOp{M, eltype(A)}(A, n)
MatrixOp{M <: AbstractMatrix}(A::M, T::Type, n::Integer) = MatrixOp{M, T}(A, n)

# MatrixOp(A::AbstractMatrix, dim_in::Tuple) = MatrixOp{typeof(A),length(dim_in)}(eltype(A), dim_in, A)
#
# MatrixOp(A::AbstractMatrix, dim_in::Vararg{Int} ) = MatrixOp(A,dim_in)
#
# MatrixOp{T<:AbstractMatrix,N}(domainType::Type, A::T, dim_in::Vararg{Int,N} ) =
# MatrixOp{T,N}(domainType,dim_in,A)
#
# MatrixOp(x::AbstractArray, A::AbstractMatrix) = MatrixOp{typeof(A),ndims(x)}(eltype(x),size(x),A)

# Mappings

A_mul_B!{M, T}(y::AbstractArray, L::MatrixOp{M, T}, b::AbstractArray) = A_mul_B!(y, L.A, b)
Ac_mul_B!{M, T}(y::AbstractArray, L::MatrixOp{M, T}, b::AbstractArray) = Ac_mul_B!(y, L.A, b)

# Properties

domainType{M, T}(L::MatrixOp{M, T}) = T
codomainType{M, T}(L::MatrixOp{M, T}) = T

function size(L::MatrixOp)
	if L.n_columns_input == 1
		( (size(L.A, 1),), (size(L.A, 2),) )
	else
		( (size(L.A, 1), L.n_columns_input), (size(L.A, 2), L.n_columns_input) )
	end
end

fun_name(L::MatrixOp) = "Matrix operator"

is_diagonal(L::MatrixOp) = isdiag(L.A)
is_full_row_rank(L::MatrixOp) = rank(L.A) == size(L.A, 1)
is_full_column_rank(L::MatrixOp) = rank(L.A) == size(L.A, 2)
# the following is O(n^3): I would assume for now no matrix is Gram diagonal
# is_gram_diagonal(L::MatrixOp)   = isdiag(L.A'*L.A)
