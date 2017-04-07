export DFT, IDFT

immutable DFT <: LinearOperator
	domainType::Type
	codomainType::Type
	dim::Tuple
	A::Base.DFT.Plan
end

size(L::DFT) = (L.dim,L.dim)

  domainType(L::DFT) = L.domainType
codomainType(L::DFT) = L.codomainType

transpose(L::DFT) = IDFT(L.codomain, L.domain, L.dim,  plan_ifft( zeros(codomainType(L),L.dim) ) )

DFT(x::AbstractArray)

immutable IDFT <: LinearOperator
	domainType::Type
	codomainType::Type
	dim::Tuple
	A::Base.DFT.Plan
end

size(L::IDFT) = (L.dim,L.dim)

  domainType(L::IDFT) = L.domainType
codomainType(L::IDFT) = L.codomainType

transpose(L::IDFT) = DFT(L.codomain, L.domain, L.dim,  plan_fft( zeros(codomainType(L),L.dim) ) )

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

