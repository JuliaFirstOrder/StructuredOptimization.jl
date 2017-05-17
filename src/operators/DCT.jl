export DCT, IDCT
abstract type CosineTransform{N,C,T1,T2} <: LinearOperator end

immutable DCT{N,
	      C<:RealOrComplex,
	      T1<:Base.DFT.Plan,
	      T2<:Base.DFT.Plan} <: CosineTransform{N,C,T1,T2}	
	dim_in::NTuple{N,Int}
	A::T1
	At::T2
end

immutable IDCT{N,
	       C<:RealOrComplex,
	       T1<:Base.DFT.Plan,
	       T2<:Base.DFT.Plan} <: CosineTransform{N,C,T1,T2}	
	dim_in::NTuple{N,Int}
	A::T1
	At::T2
end

size(L::CosineTransform) = (L.dim_in,L.dim_in)

# Constructors
function DCT{N,C<:RealOrComplex}(x::AbstractArray{C,N})  
	A,At = plan_dct(x), plan_idct(x)
	DCT{N,C,typeof(A),typeof(At)}(size(x),A,At) 
end

DCT(dim_in::Tuple) = DCT(zeros(dim_in))
DCT(T::Type,dim_in::Tuple) = DCT(zeros(T,dim_in))
DCT(dim_in::Vararg{Int64}) = DCT(dim_in)
DCT(T::Type,dim_in::Vararg{Int64}) = DCT(T,dim_in)

function IDCT{N,C<:RealOrComplex}(x::AbstractArray{C,N})  
	A,At = plan_idct(x), plan_dct(x)
	IDCT{N,C,typeof(A),typeof(At)}(size(x),A,At) 
end

IDCT(dim_in::Tuple) = IDCT(zeros(dim_in))
IDCT(T::Type,dim_in::Tuple) = IDCT(zeros(T,dim_in))
IDCT(dim_in::Vararg{Int64}) = IDCT(dim_in)
IDCT(T::Type,dim_in::Vararg{Int64}) = IDCT(T,dim_in)

# Operators

function A_mul_B!{N,C,T1,T2}(y::AbstractArray{C,N},A::DCT{N,C,T1,T2},b::AbstractArray{C,N})
	A_mul_B!(y,A.A,b)
end

function Ac_mul_B!{N,C,T1,T2}(y::AbstractArray{C,N},A::DCT{N,C,T1,T2},b::AbstractArray{C,N})
	y .= A.At*b
end

function A_mul_B!{N,C,T1,T2}(y::AbstractArray{C,N},A::IDCT{N,C,T1,T2},b::AbstractArray{C,N})
	y .= A.A*b
end

function Ac_mul_B!{N,C,T1,T2}(y::AbstractArray{C,N},A::IDCT{N,C,T1,T2},b::AbstractArray{C,N})
	A_mul_B!(y,A.At,b)
end

# Transformations


# Properties

fun_name(A::DCT)  = "Discrete Cosine Transform"
fun_name(A::IDCT) = "Inverse Discrete Cosine Transform"

domainType{N,C,T1,T2}(L::CosineTransform{N,C,T1,T2}) = C
codomainType{N,C,T1,T2}(L::CosineTransform{N,C,T1,T2}) = C

isGramDiagonal(L::CosineTransform) = true
