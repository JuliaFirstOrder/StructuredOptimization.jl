import Base: zeros

immutable Zeros{D1,D2} <: LinearOperator{D1,D2}
	x::OptVar
end
size(A::Zeros) = (size(A.x),size(A.x))

zeros{D1}(x::OptVar{D1}) = Zeros{D1,D1}(x) 

transpose{D1}(A::Zeros{D1,D1} ) = A

function A_mul_B!(y::AbstractArray,A::Zeros,b::AbstractArray)
	y .= 0
end

fun_name(A::Zeros)  = "Zeros Operator"

zeros(B::LinearOperator, args...) = NestedLinearOperator(zeros, B, args...)

immutable Empty{D1,D2} <: LinearOperator{D1,D2}
	x::OptVar
end
size(A::Empty) = (size(A.x),size(A.x))

transpose{D1,D2}(A::Empty{D1,D2}) = Empty{D2,D1}(A.x)

function A_mul_B!(y::AbstractArray,A::Empty,b::AbstractArray)
end

fun_name(A::Empty)  = ""

