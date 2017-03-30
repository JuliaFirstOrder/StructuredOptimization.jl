export diagop, DiagOp

immutable DiagOp{D1,D2} <: DiagonalOperator{D1,D2}
	sign::Bool
	d::Union{AbstractArray{D1}, Number}
	dim::Tuple

	DiagOp(d,dim) = new(true,d,dim)
	DiagOp(sign,d,dim) = new(sign,d,dim)
end
size(A::DiagOp) = (A.dim,A.dim)
-{D1,D2}(A::DiagOp{D1,D2}) = DiagOp{D1,D2}(false == sign(A),A.d,A.dim) 

DiagOp{D2}(T::Type, d::AbstractArray{D2}, dim::Vararg{Int64}) = DiagOp{T,D2}(d, dim)

DiagOp(T::Type, d::Real , dim::Vararg{Int64}) = DiagOp{T,Float64}(d, dim)
DiagOp(T::Type, d::Complex, dim::Vararg{Int64}) = DiagOp{T,Complex{Float64}}(d, dim)
DiagOp{T<:Number}(d::AbstractArray{T}, dim::Vararg{Int64}) = DiagOp(T, d, dim...)
DiagOp{T<:Number}(d::T, dim::Vararg{Int64}) = DiagOp(T, d, dim...)

function diagop{D1}(x::Variable{D1}, d) 
	A = DiagOp(D1,d, size(x)...)
	At = A'
	Affine([x], A, A', Nullable{AbstractArray}() )
end

function uA_mul_B!(y::AbstractArray,A::DiagOp,b::AbstractArray)
	y .= (*).(A.d,b)
end

function uA_mul_B!{D1<:Complex{Float64},D2<:Float64}(y::AbstractArray,A::DiagOp{D1,D2},b::AbstractArray)
	y .= real.((*).(A.d,b)) 
end

transpose{D1,D2}(A::DiagOp{D1,D2}) = DiagOp{D2,D1}(sign(A), conj(A.d),   A.dim)
inv{D1,D2}(A::DiagOp{D1,D2})       = DiagOp{D2,D1}(sign(A), (A.d).^(-1), A.dim)

fun_name(A::DiagOp)  = "Diagonal Operator"

diagop(B::AffineOperator, args...) = NestedLinearOperator(diagop,B, args...)

#other constructors with * and .*

*(d::Number, B::AffineOperator) = NestedLinearOperator(diagop, B, d)
*(d::Number, x::Variable)  = diagop(x,d)
.*(d::AbstractArray, x::Variable) = diagop(x,d)

isInvertable(A::DiagOp) = true

