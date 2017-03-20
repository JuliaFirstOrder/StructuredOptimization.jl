import Base: zeros
export Zeros

immutable Zeros{D1,D2} <: LinearOperator{D1,D2}
	dim::Tuple
end
size(A::Zeros) = (A.dim[1],A.dim[2])

Zeros(dim1::Tuple, dim2::Tuple) = Zeros{Float64,Float64}((dim1,dim2))
Zeros(T::Type, dim1::Tuple, dim2::Tuple) = Zeros{T,T}((dim1,dim2))

zeros{D1}(x::OptVar{D1}) = Affine([x], Zeros{D1,D1}((size(x),size(x))), 
				 Zeros{D1,D1}((size(x),size(x))),
				 Nullable{Vector{AbstractArray}}() )

zeros{D1}(x::OptVar{D1}, dim::Vararg{Int64}) = Affine([x], Zeros{D1,D1}((size(x),dim)), 
						     Zeros{D1,D1}((dim,size(x))),
						     Nullable{Vector{AbstractArray}}() )

transpose{D1}(A::Zeros{D1,D1} ) = Zeros{D1,D1}((A.dim[2],A.dim[1]))

function A_mul_B!(y::AbstractArray,A::Zeros,b::AbstractArray)
	y .= 0
end

fun_name(A::Zeros)  = "Zero Operator"

zeros(B::AffineOperator, args...) = NestedLinearOperator(zeros, B, args...)


