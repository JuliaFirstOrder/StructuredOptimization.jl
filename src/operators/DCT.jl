import Base: dct, idct

immutable DCT{D1,D2} <: LinearOperator{D1,D2}
	x::OptVar
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
end
size(A::DCT) = A.dim

dct{D1}(x::OptVar{D1}) = DCT{D1,D1}(x,plan_dct(x.x),plan_idct(x.x),(size(x),size(x)))

immutable IDCT{D1,D2} <: LinearOperator{D1,D2}
	x::OptVar
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
end
size(A::IDCT) = A.dim

idct{D1}(x::OptVar{D1}) = IDCT{D1,D1}(x,plan_idct(x.x),plan_dct(x.x),(size(x),size(x)))

*(A::DCT, b::AbstractArray)  = A.A*b
*(A::IDCT,b::AbstractArray)  = A.A*b

function A_mul_B!(y::AbstractArray,A::DCT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
end

function A_mul_B!(y::AbstractArray,A::IDCT,b::AbstractArray)
	#A_mul_B!(y,A.A,b) #not working??! possible bug?
	y .= A.A*b 
end

transpose{D1}(A::DCT{D1,D1} ) = IDCT{D1,D1}(A.x, A.Ainv, A.A, A.dim)
transpose{D1}(A::IDCT{D1,D1}) =  DCT{D1,D1}(A.x, A.Ainv, A.A, A.dim)

fun_name(A::DCT)  = "Discrete Cosine Transform"
fun_name(A::IDCT) = "Inverse Discrete Cosine Transform"

==(A::DCT,  B::DCT ) = A.x == B.x
==(A::IDCT, B::IDCT) = A.x == B.x

#nested Operations
dct(B::LinearOperator) = NestedLinearOperator(dct,B)
idct(B::LinearOperator) = NestedLinearOperator(idct,B)



