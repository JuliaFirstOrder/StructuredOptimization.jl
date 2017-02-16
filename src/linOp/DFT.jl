import Base: fft, ifft 

immutable DFT{D1,D2} <: LinearOp{D1,D2}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
end
size(A::DFT) = A.dim

fft{D1<:Complex}(x::OptVar{D1}) = DFT{D1,D1}(plan_fft(x.x), plan_ifft(x.x),(size(x),size(x))) 
fft{D1<:Float64}(x::OptVar{D1}) = 
DFT{D1,Complex{Float64}}(plan_fft(x.x), plan_ifft(complex(x.x)),(size(x),size(x))) 

immutable IDFT{D1,D2} <: LinearOp{D1,D2}
	A::Base.DFT.Plan
	Ainv::Base.DFT.Plan
	dim::Tuple
end
size(A::IDFT) = A.dim

ifft{D1<:Complex}(x::OptVar{D1}) = IDFT{D1,D1}(plan_ifft(x.x), plan_fft(x.x),(size(x),size(x))) 
ifft{D1<:Float64}(x::OptVar{D1}) = 
IDFT{D1,Complex{Float64}}(plan_ifft(x.x), plan_fft(complex(x.x)),(size(x),size(x))) 

*(A::DFT,b::AbstractArray)  = (A.A*b)/sqrt(length(b))
*(A::IDFT,b::AbstractArray) = (A.A*b)*sqrt(length(b))

*{D1<:Complex,D2<:Real}(A::DFT{D1,D2},b::AbstractArray)  = real((A.A*b)/sqrt(length(b)))
*{D1<:Complex,D2<:Real}(A::IDFT{D1,D2},b::AbstractArray) = real((A.A*b)*sqrt(length(b)))

function A_mul_B!(y::AbstractArray,A::DFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y ./= sqrt(length(b)) 
end

function A_mul_B!{D1<:Real,D2<:Complex}(y::AbstractArray,A::DFT{D1,D2},b::AbstractArray{D1})
	A_mul_B!(y,A.A,complex(b))
	y ./= sqrt(length(b)) 
end

function A_mul_B!{D1<:Complex,D2<:Real}(y::AbstractArray,A::DFT{D1,D2},b::AbstractArray{D1})
	y2 = complex(y) 
	A_mul_B!(y2,A.A,b)
	y .= real( (/).(y2,sqrt(length(b))) ) 
end

function A_mul_B!(y::AbstractArray,A::IDFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y .*= sqrt(length(b)) 
end

function A_mul_B!{D1<:Real,D2<:Complex}(y::AbstractArray,A::IDFT{D1,D2},b::AbstractArray{D1})
	A_mul_B!(y,A.A,complex(b))
	y .*= sqrt(length(b)) 
end

function A_mul_B!{D1<:Complex,D2<:Real}(y::AbstractArray,A::IDFT{D1,D2},b::AbstractArray{D1})
	y2 = complex(y) 
	A_mul_B!(y2,A.A,b)
	y .= real( (*).(y2,sqrt(length(b))) ) 
end

transpose{D1,D2}( A::DFT{D1,D2})  = IDFT{D2,D1}(A.Ainv, A.A, A.dim )
transpose{D1,D2}(A::IDFT{D1,D2})  =  DFT{D2,D1}(A.Ainv, A.A, A.dim )

fun_name(A::DFT) = "Discrete Fourier Transform"
fun_name(A::IDFT) = "Inverse Discrete Fourier Transform"

#nested Operations
fft(B::LinearOp) = NestedLinearOp(fft,B)
ifft(B::LinearOp) = NestedLinearOp(ifft,B)


