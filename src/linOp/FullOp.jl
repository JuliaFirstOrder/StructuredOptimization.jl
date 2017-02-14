immutable FullOp{D1} <: LinearOp{D1}
	x::OptVar{D1}
	y::OptVar{D1}
	A::AbstractArray
end

fun_name(A::FullOp)  = "Full Operator"
function *(A::AbstractMatrix, x::OptVar)  
	if size(A,1) == size(A,2)
		X = OptVar(similar(x.x))
		return FullOp(X,X,A)
	else
		return FullOp(OptVar(similar(x.x)),OptVar(zeros(typeof(x.x[1]),size(A,1))),A)
	end
end

*(A::FullOp, b::AbstractArray)  = A.A*b
transpose(A::FullOp) = FullOp(A.y,A.x,A.A')

A_mul_B!(y::AbstractArray,A::FullOp,b::AbstractArray) = A_mul_B!(y,A.A,b)  
Ac_mul_B!(y::AbstractArray,A::FullOp,b::AbstractArray) = Ac_mul_B!(y,A.A,b)  

#nested Operations
*(A::AbstractMatrix,B::LinearOp) = NestedLinearOp(*,B)
