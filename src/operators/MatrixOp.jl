export MatrixOp

immutable MatrixOp{D1,D2} <: LinearOperator{D1,D2}
	A::AbstractMatrix
end
size(A::MatrixOp) = ((size(A.A,2),),(size(A.A,1),))

MatrixOp{D2}(T::Type, A::AbstractMatrix{D2}) = MatrixOp{T,D2}(A)
MatrixOp(A::AbstractMatrix) = MatrixOp(Float64, A)

fun_name(A::MatrixOp)  = "Matrix Operator"
*{D1,D2}(A::AbstractMatrix{D2}, x::OptVar{D1}) = Affine([x], MatrixOp(D1,A), MatrixOp(D2,A'),
							Nullable{Vector{AbstractArray}}() )

*(A::MatrixOp, b::AbstractArray)  = A.A*b
transpose{D1,D2}(A::MatrixOp{D1,D2}) = MatrixOp{D2,D1}(A.A')

A_mul_B!(y::AbstractArray,A::MatrixOp,b::AbstractArray) = A_mul_B!(y,A.A,b)  
Ac_mul_B!(y::AbstractArray,A::MatrixOp,b::AbstractArray) = Ac_mul_B!(y,A.A,b)  

#nested Operations
*(A::AbstractMatrix,B::AffineOperator) = NestedLinearOperator(*, B, A)
