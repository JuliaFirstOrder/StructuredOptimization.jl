import Base: fft, ifft 

immutable DFT{T<:Union{Real,Complex}} <: LinearOp
	x::OptVar{T}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
	dom::Tuple
end
fft{T<:Real}(x::OptVar{T}) = DFT(x,plan_fft(x.x),plan_ifft(x.x),(size(x.x),size(x.x)),
										            (typeof(x.x[1]),Complex{typeof(x.x[1])})) 
fft{T<:Complex}(x::OptVar{T}) = DFT(x,plan_fft(x.x),plan_ifft(x.x),(size(x.x),size(x.x)),
										            (typeof(x.x[1]),typeof(x.x[1]))) 

immutable IDFT{T<:Union{Real,Complex}} <: LinearOp
	x::OptVar{T}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
	dom::Tuple
end
ifft{T<:Real}(x::OptVar{T}) = IDFT(x,plan_ifft(x.x),plan_fft(x.x),(size(x.x),size(x.x)),
										            (typeof(x.x[1]),Complex{typeof(x.x[1])})) 
ifft{T<:Complex}(x::OptVar{T}) = IDFT(x,plan_ifft(x.x),plan_fft(x.x),(size(x.x),size(x.x)),
										            (typeof(x.x[1]),typeof(x.x[1]))) 

*(A::DFT,b::AbstractArray)  = (A.A*b)/sqrt(length(b))
*(A::IDFT,b::AbstractArray) = (A.A*b)*sqrt(length(b))

function A_mul_B!(y::AbstractArray,A::DFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y ./= sqrt(length(b)) 
end

function A_mul_B!{T<:Real}(y::AbstractArray,A::DFT{T},b::AbstractArray)
	A_mul_B!(y,A.A,complex(b))
	y ./= sqrt(length(b)) 
end

function A_mul_B!(y::AbstractArray,A::IDFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y .*= sqrt(length(b)) 
end

function A_mul_B!{T<:Real}(y::AbstractArray,A::IDFT{T},b::AbstractArray)
	A_mul_B!(y,A.A,complex(b))
	y .*= sqrt(length(b)) 
end

transpose(A::DFT)  = IDFT(A.x, A.Ainv, A.A, A.dim, (A.dom[2],A.dom[1]))
transpose(A::IDFT) =  DFT(A.x, A.Ainv, A.A, A.dim, (A.dom[2],A.dom[1]))

fun_name(A::DFT) = "Discrete Fourier Transform"
fun_name(A::IDFT) = "Inverse Discrete Fourier Transform"
#
#
#export DFT,IDFT
#

