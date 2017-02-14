import Base: fft, ifft 

immutable DFT{D1,D2} <: LinearOp{D1,D2}
	x::OptVar{D1}
	y::OptVar{D2}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
end

function fft(x::OptVar) 
	X = OptVar(similar(x.x))
	Y = OptVar(complex(similar(x.x)))
	DFT(X, Y, plan_fft(X.x), plan_ifft(X.x)) 
end

immutable IDFT{D1,D2} <: LinearOp{D1,D2}
	x::OptVar{D1}
	y::OptVar{D2}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
end

function ifft(x::OptVar) 
	X = OptVar(similar(x.x))
	Y = OptVar(complex(similar(x.x)))
	IDFT(X, Y, plan_ifft(X.x), plan_fft(X.x)) 
end

*(A::DFT,b::AbstractArray)  = (A.A*b)/sqrt(length(b))
*(A::IDFT,b::AbstractArray) = (A.A*b)*sqrt(length(b))

*{D1<:Union{Real,Complex},D2<:Real}(A::DFT{D1,D2} ,b::AbstractArray) = real((A.A*b)/sqrt(length(b)))
*{D1<:Union{Real,Complex},D2<:Real}(A::IDFT{D1,D2},b::AbstractArray) = real((A.A*b)*sqrt(length(b)))

function A_mul_B!(y::AbstractArray,A::DFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y ./= sqrt(length(b)) 
end

function A_mul_B!{D1<:Real, D2<:Union{Real,Complex}}(y::AbstractArray,A::DFT{D1,D2},b::AbstractArray)
	A_mul_B!(y,A.A,complex(b))
	y ./= sqrt(length(b)) 
end

function A_mul_B!{D1<:Union{Real,Complex},D2<:Real}(y::AbstractArray,A::DFT{D1,D2},b::AbstractArray)
	y2 = complex(y)
	A_mul_B!(y2,A.A,b)
	y .= real(y2)./sqrt(length(b)) 
end

function A_mul_B!(y::AbstractArray,A::IDFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y .*= sqrt(length(b)) 
end

function A_mul_B!{D1<:Real, D2<:Union{Real,Complex}}(y::AbstractArray,A::IDFT{D1,D2},b::AbstractArray)
	A_mul_B!(y,A.A,complex(b))
	y .*= sqrt(length(b)) 
end

function A_mul_B!{D1<:Union{Real,Complex},D2<:Real}(y::AbstractArray,A::IDFT{D1,D2},b::AbstractArray)
	y2 = complex(y)
	A_mul_B!(y2,A.A,b)
	y .= real(y2).*sqrt(length(b)) 
end

transpose(A::DFT)  = IDFT(A.y, A.x, A.Ainv, A.A)
transpose(A::IDFT) =  DFT(A.y, A.x, A.Ainv, A.A)

fun_name(A::DFT) = "Discrete Fourier Transform"
fun_name(A::IDFT) = "Inverse Discrete Fourier Transform"

#nested Operations
fft(B::LinearOp) = NestedLinearOp(fft,B)
ifft(B::LinearOp) = NestedLinearOp(ifft,B)


