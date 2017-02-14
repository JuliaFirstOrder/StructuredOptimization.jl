import Base: reshape

immutable Reshape{D1} <: LinearOp{D1}
	x::OptVar{D1}
	y::OptVar{D1}
end

function reshape(x::OptVar, args...) 
	X = OptVar(similar(x.x))
	Reshape(X, OptVar(reshape(X.x,args...)  ))
end

*(A::Reshape,b::AbstractArray)  = reshape(b,size(A)[2])

function A_mul_B!(y::AbstractArray,A::Reshape,b::AbstractArray) 
	copy!(y, reshape(b, size(A)[2]))
end

transpose(A::Reshape) = Reshape(A.y,A.x)

fun_name(A::Reshape) = "Reshape Operator"

#nested Operations
reshape(B::LinearOp,args...) = NestedLinearOp(reshape,B, args...)
