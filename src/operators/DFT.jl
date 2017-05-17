export DFT, IDFT

abstract type FourierTransform{N,C,D,T1,T2} <: LinearOperator end

immutable DFT{N,
	      C<:RealOrComplex,
	      D<:RealOrComplex,
	      T1<:Base.DFT.Plan,
	      T2<:Base.DFT.Plan} <: FourierTransform{N,C,D,T1,T2}
	dim_in::NTuple{N,Int}
	A::T1
	At::T2
end

immutable IDFT{N,
	       C<:RealOrComplex,
	       D<:RealOrComplex,
	       T1<:Base.DFT.Plan,
	       T2<:Base.DFT.Plan} <: FourierTransform{N,C,D,T1,T2}
	dim_in::NTuple{N,Int}
	A::T1
	At::T2
end

# Constructors

function DFT{N,D<:Real}(x::AbstractArray{D,N})
	A,At = plan_fft(x), plan_bfft(fft(x))
	DFT{N,Complex{D},D,typeof(A),typeof(At)}(size(x),A,At)
end

function DFT{N,D<:Complex}(x::AbstractArray{D,N})
	A,At = plan_fft(x), plan_bfft(fft(x))
	DFT{N,D,D,typeof(A),typeof(At)}(size(x),A,At)
end

DFT(dim_in::Tuple) = DFT(zeros(dim_in))
DFT(T::Type,dim_in::Tuple) = DFT(zeros(T,dim_in))
DFT(dim_in::Vararg{Int}) = DFT(dim_in)
DFT(T::Type,dim_in::Vararg{Int}) = DFT(T,dim_in)

function IDFT{N,D<:Real}(x::AbstractArray{D,N})
	A,At = plan_ifft(x), plan_fft(ifft(x))
	IDFT{N,Complex{D},D,typeof(A),typeof(At)}(size(x),A,At)
end

function IDFT{N,D<:Complex}(x::AbstractArray{D,N})
	A,At = plan_ifft(x), plan_fft(ifft(x))
	IDFT{N,D,D,typeof(A),typeof(At)}(size(x),A,At)
end

IDFT(dim_in::Tuple) = IDFT(zeros(dim_in))
IDFT(T::Type,dim_in::Tuple) = IDFT(zeros(T,dim_in))
IDFT(dim_in::Vararg{Int}) = IDFT(dim_in)
IDFT(T::Type,dim_in::Vararg{Int}) = IDFT(T,dim_in)

# Mappings

function A_mul_B!{N,C<:Complex,T1,T2}(y::AbstractArray{C,N},
				      L::DFT{N,C,C,T1,T2},
				      b::AbstractArray{C,N})
	A_mul_B!(y,L.A,b)
end

function A_mul_B!{N,C<:Complex,D<:Real,T1,T2}(y::AbstractArray{C,N},
					      L::DFT{N,C,D,T1,T2},
					      b::AbstractArray{D,N})
	A_mul_B!(y,L.A,complex(b))
end

function Ac_mul_B!{N,C<:Complex,T1,T2}(y::AbstractArray{C,N},
				       L::DFT{N,C,C,T1,T2},
				       b::AbstractArray{C,N})
	A_mul_B!(y,L.At,b)
end

function Ac_mul_B!{N,C<:Complex,D<:Real,T1,T2}(y::AbstractArray{D,N},
					      L::DFT{N,C,D,T1,T2},
					      b::AbstractArray{C,N})
	y2 = complex(y)
	A_mul_B!(y2,L.At,b)
	y .= real.(y2)
end

function A_mul_B!{N,C<:Complex,T1,T2}(y::AbstractArray{C,N},
				      L::IDFT{N,C,C,T1,T2},
				      b::AbstractArray{C,N})
	A_mul_B!(y,L.A,b)
end

function A_mul_B!{N,C<:Complex,D<:Real,T1,T2}(y::AbstractArray{C,N},
					      L::IDFT{N,C,D,T1,T2},
					      b::AbstractArray{D,N})
	A_mul_B!(y,L.A,complex(b))
end

function Ac_mul_B!{N,C<:Complex,T1,T2}(y::AbstractArray{C,N},
				       L::IDFT{N,C,C,T1,T2},
				       b::AbstractArray{C,N})
	A_mul_B!(y,L.At,b)
	y ./= length(b)
end

function Ac_mul_B!{N,C<:Complex,D<:Real,T1,T2}(y::AbstractArray{D,N},
					       L::IDFT{N,C,D,T1,T2},
					       b::AbstractArray{C,N})

	y2 = complex(y)
	A_mul_B!(y2,L.At,b)
	y .= (/).(real.(y2), length(b))
end

# Properties

size(L::FourierTransform) = (L.dim_in,L.dim_in)

fun_name(A::DFT) = "Discrete Fourier Transform"
fun_name(A::IDFT) = "Inverse Discrete Fourier Transform"

domainType{N,C,D,T1,T2}(L::FourierTransform{N,C,D,T1,T2}) = D
codomainType{N,C,D,T1,T2}(L::FourierTransform{N,C,D,T1,T2}) = C

is_gram_diagonal(L::FourierTransform) = true
