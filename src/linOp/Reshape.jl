import Base: reshape

immutable Reshape{T} <: LinearOp
	x::OptVar{T}
	dim::Tuple
	dom::Tuple
end

reshape(x::OptVar, newdim::Vararg{Int64}) = Reshape(x, (size(x),newdim), 
																									   	 (typeof(x.x[1]),typeof(x.x[1])))
reshape(x::OptVar, newdim::Tuple) = Reshape(x, (size(x),newdim), 
																									   	 (typeof(x.x[1]),typeof(x.x[1])))
*(A::Reshape,b::AbstractArray)  = reshape(b,A.dim[2])

function A_mul_B!(y::AbstractArray,A::Reshape,b::AbstractArray) 
	copy!(y, reshape(b, A.dim[2]))
end

transpose(A::Reshape) = Reshape(OptVar(reshape(A.x.x,A.dim[2])),(A.dim[2],A.dim[1]),A.dom)

export Reshape

fun_name(A::Reshape) = "Reshape Operator"

#nested Operations
reshape(B::LinearOp,args...) = NestedLinearOp(reshape,B, args...)
