import Base: dct, idct

immutable DCT{T} <: LinearOp
	x::OptVar{T}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
	dom::Tuple
end
dct(x::OptVar) = DCT(x,plan_dct(x.x),
										 plan_idct(x.x),(size(x.x),size(x.x)),(typeof(x.x[1]),typeof(x.x[1])))

immutable IDCT{T} <: LinearOp
	x::OptVar{T}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
	dom::Tuple
end
idct(x::OptVar) = IDCT(x,plan_idct(x.x),
										   plan_dct(x.x),(size(x.x),size(x.x)),(typeof(x.x[1]),typeof(x.x[1])))

*(A::DCT,b::AbstractArray)  = A.A*b
*(A::IDCT,b::AbstractArray) = A.A*b

function A_mul_B!(y::AbstractArray,A::DCT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
end

function A_mul_B!(y::AbstractArray,A::IDCT,b::AbstractArray)
	#A_mul_B!(y,A.A,b) #not working??! possible bug?
	y .= A.A*b 
end

transpose(A::DCT)  = IDCT(A.x, A.Ainv, A.A, A.dim, A.dom )
transpose(A::IDCT) =  DCT(A.x, A.Ainv, A.A, A.dim, A.dom )

fun_name(A::DCT)  = "Discrete Cosine Transform"
fun_name(A::IDCT) = "Inverse Cosine Fourier Transform"

#nested Operations
dct(B::LinearOp) = NestedLinearOp(dct,B)
idct(B::LinearOp) = NestedLinearOp(idct,B)



