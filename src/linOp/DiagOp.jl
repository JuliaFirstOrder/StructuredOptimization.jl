
immutable DiagOp{D1,D2} <: LinearOp{D1,D2}
	x::OptVar
	d::Union{AbstractArray{D1}, Number}
	dim::Tuple
end
size(A::DiagOp) = A.dim

diagop{D1}(x::OptVar{D1}, d::AbstractArray{D1}) = DiagOp{D1,D1}(x,d,(size(x),size(x)))
diagop{D1,T<:Number}(x::OptVar{D1}, d::T) = DiagOp{D1,D1}(x,d,(size(x),size(x)))

function A_mul_B!{T}(y::AbstractArray{T},A::DiagOp,b::AbstractArray{T})
	y .= (*).(A.d,y)
	copy!(y,1,b,1,length(b))
end

transpose{D1}(A::DiagOp{D1,D1}) = DiagOp{D1,D1}(A.x,conj(A.d),A.dim)
fun_name(A::DiagOp)  = "Diagonal Operator"

diagop(B::LinearOp, args...) = NestedLinearOp(diagop,B, args...)

export diagop
