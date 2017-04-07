export MatrixOp

immutable MatrixOp <: LinearOperator
	A::AbstractMatrix
	size_in::Tuple
end

MatrixOp(A::AbstractMatrix) = MatrixOp(A, (size(A, 2),))
MatrixOp(A::AbstractMatrix, x::AbstractArray) = MatrixOp(A, size(x))

size(L::MatrixOp) = length(L.size_in) == 1 ? ((size(L.A,1),), L.size_in) : ((size(L.A,1), L.size_in[2]), L.size_in)

A_mul_B!(y, A::MatrixOp, b) = A_mul_B!(y, A.A, b)
At_mul_B!(y, A::MatrixOp, b) = At_mul_B!(y, A.A, b)

fun_name(A::MatrixOp)  = "Matrix operator"

# function *(A::MatrixOp,b::AbstractArray)
# 	C = codomainType(A)
# 	y = Array{C}(size(A,1))
# 	A_mul_B!(y,A,b)
# 	return y
# end

################################################################################
# FROM HERE ON IT IS USERS' SYNTAX
################################################################################

*{D1,D2}(A::AbstractMatrix{D2}, x::Variable{D1}) = Affine([x], MatrixOp(D1,A), MatrixOp(D2,A'), Nullable{AbstractArray}() )

#nested Operations
*(A::AbstractMatrix,B::AffineOperator) = NestedLinearOperator(*, B, A)
