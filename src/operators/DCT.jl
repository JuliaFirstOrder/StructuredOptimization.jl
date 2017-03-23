import Base: dct, idct
export DCT, IDCT

immutable DCT{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	A::Base.DFT.Plan
	dim::Tuple

	DCT(A,dim)      = new(true,A,dim)
	DCT(sign,A,dim) = new(sign,A,dim)
end
size(A::DCT) = (A.dim, A.dim)
-{D1,D2}(A::DCT{D1,D2}) = DCT{D1,D2}(false == sign(A), A.A, A.dim) 

DCT{D1}(x::AbstractArray{D1}) = DCT{D1,D1}(plan_dct(x), size(x))
transpose{D1}(A::DCT{D1,D1})  = IDCT{D1,D1}(sign(A), plan_idct( zeros(D1,A.dim) ), A.dim )

function dct(x::OptVar)  
	A = DCT(x.x) 
	Affine([x], A, A', Nullable{AbstractArray}() )
end

immutable IDCT{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	A::Base.DFT.Plan
	dim::Tuple

	IDCT(A,dim)      = new(true,A,dim)
	IDCT(sign,A,dim) = new(sign,A,dim)
end
size(A::IDCT) = (A.dim, A.dim)
-{D1,D2}(A::IDCT{D1,D2}) = IDCT{D1,D2}(false == sign(A), A.A, A.dim) 

IDCT{D1}(x::AbstractArray{D1}) = IDCT{D1,D1}(plan_idct(x), size(x))
transpose{D1}(A::IDCT{D1,D1})  =  DCT{D1,D1}(sign(A), plan_dct( zeros(D1,A.dim) ), A.dim )

function idct(x::OptVar)  
	A = IDCT(x.x) 
	Affine([x], A, A', Nullable{AbstractArray}() )
end

function uA_mul_B!(y::AbstractArray,A::DCT,b::AbstractArray)
	A_mul_B!(y,A.A,b) 
end

function uA_mul_B!(y::AbstractArray,A::IDCT,b::AbstractArray)
	#A_mul_B!(y,A.A,b) #not working??! possible bug?
	y .= A.A*b 
end

fun_name(A::DCT)  = "Discrete Cosine Transform"
fun_name(A::IDCT) = "Inverse Discrete Cosine Transform"

#nested Operations
dct(B::AffineOperator)  = NestedLinearOperator(dct,B)
idct(B::AffineOperator) = NestedLinearOperator(idct,B)



