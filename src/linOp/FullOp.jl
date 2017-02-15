immutable FullOp{D1,D2} <: LinearOp{D1,D2}
	A::AbstractArray
	dim::Tuple
end

fun_name(A::FullOp)  = "Full Operator"
*{D1}(A::AbstractMatrix, x::OptVar{D1}) = FullOp{D1,D1}(A, ((size(A,2),),(size(A,1),)) ) 

*(A::FullOp, b::AbstractArray)  = A.A*b
transpose{D1}(A::FullOp{D1,D1}) = FullOp{D1,D1}(A.A',(A.dim[2],A.dim[1]))

A_mul_B!(y::AbstractArray,A::FullOp,b::AbstractArray) = A_mul_B!(y,A.A,b)  
Ac_mul_B!(y::AbstractArray,A::FullOp,b::AbstractArray) = Ac_mul_B!(y,A.A,b)  

#nested Operations
*(A::AbstractMatrix,B::LinearOp) = NestedLinearOp(*, B, A)
