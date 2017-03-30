export Empty, emptyop

immutable Empty{D1,D2} <: LinearOperator{D1,D2}
	dim::Tuple
end
size(A::Empty) = (A.dim[1], A.dim[2])

emptyop{D1}(x::Variable{D1}) = Affine([x], Empty{D1,D1}((size(x),size(x))), 
				   Empty{D1,D1}((size(x),size(x))),
				   Nullable{Vector{AbstractArray}}() )

#creates an empty operator from another operator
emptyop{D1,D2}(A::LinearOperator{D1,D2}) = Empty{D1,D2}(size(A))

transpose{D1,D2}(A::Empty{D1,D2}) = Empty{D2,D1}((A.dim[2],A.dim[1]))

function A_mul_B!(y::AbstractArray,A::Empty,b::AbstractArray)
end

fun_name(A::Empty)  = ""
