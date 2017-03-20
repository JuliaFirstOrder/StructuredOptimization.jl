import Base: fft, ifft 
export DFT, IDFT

immutable DFT{D1,D2} <: LinearOperator{D1,D2}
	A::Base.DFT.Plan
	dim::Tuple
end
size(A::DFT) = (A.dim, A.dim)

DFT{D1}(x::AbstractArray{D1}) = DFT{D1,Complex{Float64}}(plan_fft(x), size(x))
transpose{D1,D2}(A::DFT{D1,D2})  = IDFT{D2,D1}( plan_ifft( zeros(D2,A.dim) ), A.dim )

function fft(x::OptVar)  
	A = DFT(x.x) 
	Affine([x], A, A', Nullable{Vector{AbstractArray}}() )
end

immutable IDFT{D1,D2} <: LinearOperator{D1,D2}
	A::Base.DFT.Plan
	dim::Tuple
end
size(A::IDFT) = (A.dim, A.dim)

IDFT{D1}(x::AbstractArray{D1}) = IDFT{D1,Complex{Float64}}(plan_ifft(x), size(x))
transpose{D1,D2}(A::IDFT{D1,D2})  = DFT{D2,D1}( plan_fft( zeros(D2,A.dim) ), A.dim )

function ifft(x::OptVar)  
	A = IDFT(x.x) 
	Affine([x], A, A', Nullable{Vector{AbstractArray}}() )
end

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
	y .= real.( (/).(y2,sqrt(length(b))) ) 
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
	y .= real.( (*).(y2,sqrt(length(b))) ) 
end

fun_name(A::DFT) = "Discrete Fourier Transform"
fun_name(A::IDFT) = "Inverse Discrete Fourier Transform"

#nested Operations
fft(B::AffineOperator) = NestedLinearOperator(fft,B)
ifft(B::AffineOperator) = NestedLinearOperator(ifft,B)


