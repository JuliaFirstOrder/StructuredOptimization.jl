function compose(f::Type, B::Affine, args...)
	mid = Array{codomainType(operator(B))}(size(operator(B),1))
	A = f(mid, args...)
	C = Compose(A,operator(B),mid)
	Affine(variable(B),C)
end

function compose(f::Type, B::TiltedAffine, args...)
	mid = Array{codomainType(operator(B))}(size(operator(B),1))
	A = f(mid, args...)
	C = Compose(A,operator(B),mid)
	Affine(variable(B),C)+A*tilt(B)
end

(*)(L::LinearOperator,A::Affine)         = Affine(variable(A),L*operator(A))
(*)(L::LinearOperator,A::TiltedAffine)   = Affine(variable(A),L*operator(A))+L*tilt(A)

(*){T<:Number}(coeff::T,A::Affine)       = Affine(variable(A),coeff*operator(A))
(*){T<:Number}(coeff::T,A::TiltedAffine) = Affine(variable(A),coeff*operator(A))+coeff.*tilt(A)

(*){T<:Number}(coeff::T,x::Variable)     = Affine(x,coeff*Eye(size(x)))

import Base: zeros, 
             eye,
	     reshape,
	     .*,
	     getindex,
	     *,
	     fft,
	     ifft,
	     dct,
	     idct,
	     conv,
	     xcorr

export       finitediff,
	     variation
       
	     
zeros(x::AbstractVariable) = Zeros(~x)*x
zeros(A::AffineOperator, args...) = compose(Zeros, A)

eye(x::AbstractVariable) = Eye(~x)*x
eye(A::AffineOperator  ) = A

reshape(x::AbstractVariable, args...) = Reshape(Eye(~x), args...)*x
reshape(A::Affine          , args...) = Affine(variable(A), Reshape(operator(A), args...))
reshape(A::TiltedAffine    , args...) = 
Affine(variable(A), Reshape(operator(A), args...))+reshape(tilt(A), args...)

(.*)(d::AbstractVector,x::AbstractVariable) =   DiagOp(~x, d)*x
(.*)(d::AbstractVector,A::AffineOperator  ) =   compose(DiagOp, A, d)

getindex(x::AbstractVariable, args...) = GetIndex(~x, args)*x
getindex(A::AffineOperator,   args...) = compose(GetIndex, A, args)

(*)(M::AbstractMatrix,x::AbstractVariable)  = MatrixOp(~x, M)*x
(*)(M::AbstractMatrix,A::AffineOperator  )  = compose(MatrixOp, A, M)

fft(x::AbstractVariable) = DFT(~x)*x
fft(A::AffineOperator  ) = compose(DFT, A)

ifft(x::AbstractVariable) = IDFT(~x)*x
ifft(A::AffineOperator  ) = compose(IDFT, A)

dct(x::AbstractVariable) = DCT(~x)*x
dct(A::AffineOperator  ) = compose(DCT, A)

idct(x::AbstractVariable) = IDCT(~x)*x
idct(A::AffineOperator  ) = compose(IDCT, A)

finitediff(x::AbstractVariable, args...) = FiniteDiff(~x, args...)*x
finitediff(A::AffineOperator,   args...) = compose(FiniteDiff, A, args...)

variation(x::AbstractVariable, args...) = Variation(~x, args...)*x
variation(A::AffineOperator,   args...) = compose(Variation, A, args...)

conv(x::AbstractVariable, args...)  = Conv(~x, args...)*x
conv(A::AffineOperator,   args...) = compose(Conv, A, args...)

xcorr(x::AbstractVariable, args...) = Xcorr(~x, args...)*x
xcorr(A::AffineOperator,   args...) = compose(Xcorr, A, args...)


