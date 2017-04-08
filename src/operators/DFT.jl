export DFT, IDFT
abstract FourierTransform <: LinearOperator

immutable DFT <: FourierTransform
	dim_in::Tuple
	domainType::Type
	codomainType::Type
	A::Base.DFT.Plan
end

immutable IDFT <: FourierTransform
	dim_in::Tuple
	domainType::Type
	codomainType::Type
	A::Base.DFT.Plan
end

size(L::FourierTransform) = (L.dim_in,L.dim_in)

  domainType(L::FourierTransform) = L.domainType
codomainType(L::FourierTransform) = L.codomainType

DFT(x::AbstractArray)  =  DFT(size(x),eltype(x),Complex{Float64},plan_fft(x)) 
IDFT(x::AbstractArray) = IDFT(size(x),eltype(x),Complex{Float64},plan_ifft(x)) 

transpose(L::DFT) = IDFT(L.dim_in, L.codomainType, L.domainType, plan_ifft( zeros(codomainType(L),L.dim_in) ) )

transpose(L::IDFT) = DFT(L.dim_in, L.codomainType, L.domainType, plan_fft( zeros(codomainType(L),L.dim_in) ) )


function A_mul_B!(y::AbstractArray,L::DFT,b::AbstractArray)
	if domainType(L) == codomainType(L)
		A_mul_B!(y,L.A,b)
		y ./= sqrt(length(b))
	elseif domainType(L) <: Real && codomainType(L) <: Complex
		A_mul_B!(y,L.A,complex(b))
		y ./= sqrt(length(b))
	elseif domainType(L) <: Complex && codomainType(L) <: Real
		y2 = complex(y)
		A_mul_B!(y2,L.A,b)
		y .= real.( (/).(y2,sqrt(length(b))) )
	end
end

function A_mul_B!(y::AbstractArray,L::IDFT,b::AbstractArray)
	if domainType(L) == codomainType(L)
		A_mul_B!(y,L.A,b)
		y .*= sqrt(length(b))
	elseif domainType(L) <: Real && codomainType(L) <: Complex
		A_mul_B!(y,L.A,complex(b))
		y .*= sqrt(length(b))
	elseif domainType(L) <: Complex && codomainType(L) <: Real
		y2 = complex(y)
		A_mul_B!(y2,L.A,b)
		y .= real.( (*).(y2,sqrt(length(b))) )
	end
end

fun_name(A::DFT) = "Discrete Fourier Transform"
fun_name(A::IDFT) = "Inverse Discrete Fourier Transform"

