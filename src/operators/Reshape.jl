import Base: reshape
export Reshape

immutable Reshape{D1,D2} <: LinearOperator{D1,D2}
	dim::Tuple
end
size(A::Reshape) = (A.dim[1],A.dim[2])
Reshape{D1}(x::AbstractArray{D1}, dim...) = Reshape{D1,D1}((size(x),dim))

function reshape{D1}(x::OptVar{D1}, dim::Vararg{Int64}) 
	A = Reshape(x.x, dim...)
	Affine([x], A, A',Nullable{Vector{AbstractArray}}() )
end

*(A::Reshape,b::AbstractArray)  = reshape(b,A.dim[2]...)

function A_mul_B!(y::AbstractArray,A::Reshape,b::AbstractArray) 
	y .= reshape(b,A.dim[2]...)
end

transpose{D1}(A::Reshape{D1,D1}) = Reshape{D1,D1}((A.dim[2],A.dim[1]))

fun_name(A::Reshape) = "Reshape Operator"

#nested Operations
reshape(B::AffineOperator,args...) = NestedLinearOperator(reshape,B, args...)
