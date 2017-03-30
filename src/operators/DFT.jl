import Base: fft, ifft 
export DFT, IDFT

immutable DFT{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	A::Base.DFT.Plan
	dim::Tuple

	DFT(A,dim)      = new(true,A,dim)
	DFT(sign,A,dim) = new(sign,A,dim)
end
size(A::DFT) = (A.dim, A.dim)
-{D1,D2}(A::DFT{D1,D2}) = DFT{D1,D2}(false == sign(A), A.A, A.dim) 

DFT{D1}(x::AbstractArray{D1}) = DFT{D1,Complex{Float64}}(plan_fft(x), size(x))
transpose{D1,D2}(A::DFT{D1,D2})  = IDFT{D2,D1}(sign(A), plan_ifft( zeros(D2,A.dim) ), A.dim )

function fft(x::Variable)  
	A = DFT(x.x) 
	Affine([x], A, A', Nullable{AbstractArray}() )
end

immutable IDFT{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	A::Base.DFT.Plan
	dim::Tuple

	IDFT(A,dim)      = new(true,A,dim)
	IDFT(sign,A,dim) = new(sign,A,dim)
end
size(A::IDFT) = (A.dim, A.dim)
-{D1,D2}(A::IDFT{D1,D2}) = IDFT{D1,D2}(false == sign(A), A.A, A.dim) 

IDFT{D1}(x::AbstractArray{D1})    = IDFT{D1,Complex{Float64}}(plan_ifft(x), size(x))
transpose{D1,D2}(A::IDFT{D1,D2})  = DFT{D2,D1}(sign(A), plan_fft( zeros(D2,A.dim) ), A.dim )

function ifft(x::Variable)  
	A = IDFT(x.x) 
	Affine([x], A, A', Nullable{AbstractArray}() )
end

function uA_mul_B!(y::AbstractArray,A::DFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y ./= sqrt(length(b)) 
end

function uA_mul_B!{D1<:Real,D2<:Complex}(y::AbstractArray,A::DFT{D1,D2},b::AbstractArray{D1})
	A_mul_B!(y,A.A,complex(b))
	y ./= sqrt(length(b)) 
end

function uA_mul_B!{D1<:Complex,D2<:Real}(y::AbstractArray,A::DFT{D1,D2},b::AbstractArray{D1})
	y2 = complex(y) 
	A_mul_B!(y2,A.A,b)
	y .= real.( (/).(y2,sqrt(length(b))) ) 
end

function uA_mul_B!(y::AbstractArray,A::IDFT,b::AbstractArray)
	A_mul_B!(y,A.A,b)
	y .*= sqrt(length(b)) 
end

function uA_mul_B!{D1<:Real,D2<:Complex}(y::AbstractArray,A::IDFT{D1,D2},b::AbstractArray{D1})
	A_mul_B!(y,A.A,complex(b))
	y .*= sqrt(length(b)) 
end

function uA_mul_B!{D1<:Complex,D2<:Real}(y::AbstractArray,A::IDFT{D1,D2},b::AbstractArray{D1})
	y2 = complex(y) 
	A_mul_B!(y2,A.A,b)
	y .= real.( (*).(y2,sqrt(length(b))) ) 
end

fun_name(A::DFT) = "Discrete Fourier Transform"
fun_name(A::IDFT) = "Inverse Discrete Fourier Transform"

#nested Operations
fft(B::AffineOperator) = NestedLinearOperator(fft,B)
ifft(B::AffineOperator) = NestedLinearOperator(ifft,B)


