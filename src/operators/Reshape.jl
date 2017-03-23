import Base: reshape
export Reshape

immutable Reshape{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	dim::Tuple

	Reshape(sign,dim) = new(sign,dim)
	Reshape(dim)      = new(true,dim)
end
size(A::Reshape) = (A.dim[1],A.dim[2])
-{D1,D2}(A::Reshape{D1,D2}) = Reshape{D1,D2}(false == sign(A), A.dim) 

Reshape{D1}(x::AbstractArray{D1}, dim...) = Reshape{D1,D1}((size(x),dim))

function reshape{D1}(x::OptVar{D1}, dim::Vararg{Int64}) 
	A = Reshape(x.x, dim...)
	Affine([x], A, A',Nullable{Vector{AbstractArray}}() )
end

function uA_mul_B!(y::AbstractArray,A::Reshape,b::AbstractArray) 
	y .= reshape(b,A.dim[2]...)
end

transpose{D1}(A::Reshape{D1,D1}) = Reshape{D1,D1}(sign(A),(A.dim[2],A.dim[1]))

fun_name(A::Reshape) = "Reshape Operator"

#nested Operations
reshape(B::AffineOperator,args...) = NestedLinearOperator(reshape,B, args...)
