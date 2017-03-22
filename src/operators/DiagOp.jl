export diagop, DiagOp

immutable DiagOp{D1,D2} <: DiagonalOperator{D1,D2}
	d::Union{AbstractArray{D1}, Number}
	dim::Tuple
end
size(A::DiagOp) = (A.dim,A.dim)

DiagOp{D2}(T::Type, d::AbstractArray{D2}, dim::Vararg{Int64}) = DiagOp{T,D2}(d, dim)

DiagOp(T::Type, d::Real , dim::Vararg{Int64}) = DiagOp{T,Float64}(d, dim)
DiagOp(T::Type, d::Complex, dim::Vararg{Int64}) = DiagOp{T,Complex{Float64}}(d, dim)
DiagOp{T<:Number}(d::AbstractArray{T}, dim::Vararg{Int64}) = DiagOp(T, d, dim...)
DiagOp{T<:Number}(d::T, dim::Vararg{Int64}) = DiagOp(T, d, dim...)

function diagop{D1}(x::OptVar{D1}, d) 
	A = DiagOp(D1,d, size(x)...)
	At = A'
	Affine([x], A, A', Nullable{Vector{AbstractArray}}() )
end

function A_mul_B!(y::AbstractArray,A::DiagOp,b::AbstractArray)
	y .= (*).(A.d,b)
end

function A_mul_B!{D1<:Complex{Float64},D2<:Float64}(y::AbstractArray,A::DiagOp{D1,D2},b::AbstractArray)
	y .= real.((*).(A.d,b))
end

transpose{D1,D2}(A::DiagOp{D1,D2}) = DiagOp{D2,D1}(conj(A.d),   A.dim)
inv{D1,D2}(A::DiagOp{D1,D2})       = DiagOp{D2,D1}((A.d).^(-1), A.dim)

fun_name(A::DiagOp)  = "Diagonal Operator"

diagop(B::AffineOperator, args...) = NestedLinearOperator(diagop,B, args...)

#other constructors with * and .*

*(d::Number, B::AffineOperator) = NestedLinearOperator(diagop, B, d)
*(d::Number, x::OptVar)  = diagop(x,d)
.*(d::AbstractArray, x::OptVar) = diagop(x,d)

isInvertable(A::DiagOp) = true

