
immutable DiagOp{D1,D2} <: DiagonalOperator{D1,D2}
	x::OptVar
	d::Union{AbstractArray{D1}, Number}
end
size(A::DiagOp) = (size(A.x),size(A.x))

diagop{D1,D2}(x::OptVar{D1}, d::AbstractArray{D2}) = DiagOp{D1,D2}(x,d)
diagop{D1,T<:Number}(x::OptVar{D1}, d::T) = DiagOp{D1,D1}(x,d)
diagop{D1,T<:Complex}(x::OptVar{D1}, d::T) = DiagOp{D1,Complex{Float64}}(x,d)

function A_mul_B!(y::AbstractArray,A::DiagOp,b::AbstractArray)
	y .= (*).(A.d,b)
end

function A_mul_B!{D1<:Complex{Float64},D2<:Float64}(y::AbstractArray,A::DiagOp{D1,D2},b::AbstractArray)
	y .= real.((*).(A.d,b))
end

transpose{D1,D2}(A::DiagOp{D1,D2}) = DiagOp{D2,D1}(A.x,conj(A.d))
inv{D1,D2}(A::DiagOp{D1,D2})       = DiagOp{D2,D1}(A.x,(A.d).^(-1))

fun_name(A::DiagOp)  = "Diagonal Operator"

diagop(B::LinearOperator, args...) = NestedLinearOperator(diagop,B, args...)
*(d::Union{Float64,Complex{Float64}}, B::LinearOperator) = NestedLinearOperator(diagop, B, d)
*{D1,T<:Number}(d::T, x::OptVar{D1})  = DiagOp{D1,D1}(x,d)
*{D1,T<:Complex}(d::T, x::OptVar{D1}) = DiagOp{D1,Complex{Float64}}(x,d)

isInvertable(A::DiagOp) = true

export diagop
