immutable FullOp{D1,D2} <: LinearOperator{D1,D2}
	x::OptVar
	A::AbstractArray
end
size(A::FullOp) = ((size(A.A,2),),(size(A.A,1),))

fun_name(A::FullOp)  = "Full Operator"
*{D1}(A::AbstractMatrix, x::OptVar{D1}) = FullOp{D1,D1}(x, A) 

*(A::FullOp, b::AbstractArray)  = A.A*b
transpose{D1}(A::FullOp{D1,D1}) = FullOp{D1,D1}(A.x, A.A')

A_mul_B!(y::AbstractArray,A::FullOp,b::AbstractArray) = A_mul_B!(y,A.A,b)  
Ac_mul_B!(y::AbstractArray,A::FullOp,b::AbstractArray) = Ac_mul_B!(y,A.A,b)  

#nested Operations
*(A::AbstractMatrix,B::LinearOperator) = NestedLinearOperator(*, B, A)
