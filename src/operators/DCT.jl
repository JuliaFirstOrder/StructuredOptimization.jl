import Base: dct, idct
export DCT, IDCT

immutable DCT{D1,D2} <: LinearOperator{D1,D2}
	A::Base.DFT.Plan
	dim::Tuple
end
size(A::DCT) = (A.dim, A.dim)

DCT{D1}(x::AbstractArray{D1}) = DCT{D1,D1}(plan_dct(x), size(x))
transpose{D1}(A::DCT{D1,D1})  = IDCT{D1,D1}( plan_idct( zeros(D1,A.dim) ), A.dim )

function dct(x::OptVar)  
	A = DCT(x.x) 
	Affine([x], A, A', Nullable{Vector{AbstractArray}}() )
end

immutable IDCT{D1,D2} <: LinearOperator{D1,D2}
	A::Base.DFT.Plan
	dim::Tuple
end
size(A::IDCT) = (A.dim, A.dim)

IDCT{D1}(x::AbstractArray{D1}) = IDCT{D1,D1}(plan_idct(x), size(x))
transpose{D1}(A::IDCT{D1,D1})  =  DCT{D1,D1}( plan_dct( zeros(D1,A.dim) ), A.dim )

function idct(x::OptVar)  
	A = IDCT(x.x) 
	Affine([x], A, A', Nullable{Vector{AbstractArray}}() )
end

*(A::DCT, b::AbstractArray)  = A.A*b
*(A::IDCT,b::AbstractArray)  = A.A*b

function A_mul_B!(y::AbstractArray,A::DCT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
end

function A_mul_B!(y::AbstractArray,A::IDCT,b::AbstractArray)
	#A_mul_B!(y,A.A,b) #not working??! possible bug?
	y .= A.A*b 
end

fun_name(A::DCT)  = "Discrete Cosine Transform"
fun_name(A::IDCT) = "Inverse Discrete Cosine Transform"

#nested Operations
dct(B::AffineOperator)  = NestedLinearOperator(dct,B)
idct(B::AffineOperator) = NestedLinearOperator(idct,B)



