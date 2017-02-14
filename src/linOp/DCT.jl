import Base: dct, idct

immutable DCT{D1} <: LinearOp{D1}
	x::OptVar{D1}
	y::OptVar{D1}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
end

function dct(x::OptVar) 
	X = OptVar(similar(x.x))
	DCT(X,X,plan_dct(x.x),plan_idct(x.x))
end

immutable IDCT{D1} <: LinearOp{D1}
	x::OptVar{D1}
	y::OptVar{D1}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
end
function idct(x::OptVar) 
	X = OptVar(similar(x.x))
	IDCT(X,X,plan_idct(x.x),plan_dct(x.x))
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

transpose(A::DCT)  = IDCT(A.y, A.x, A.Ainv, A.A)
transpose(A::IDCT) =  DCT(A.y, A.x, A.Ainv, A.A)

fun_name(A::DCT)  = "Discrete Cosine Transform"
fun_name(A::IDCT) = "Inverse Discrete Cosine Transform"

#nested Operations
dct(B::LinearOp) = NestedLinearOp(dct,B)
idct(B::LinearOp) = NestedLinearOp(idct,B)



