
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


(.*)(d::AbstractVector,A::AffineOperator  ) =   compose(DiagOp, A, d)

(*)(M::AbstractMatrix,A::AffineOperator  )  = compose(MatrixOp, A, M)

eye(       A::AffineOperator         ) = A

zeros(     A::AffineOperator, args...) = compose(Zeros,      A, args...)

getindex(  A::AffineOperator, args...) = compose(GetIndex,   A, args)

fft(       A::AffineOperator, args...) = compose(DFT,        A, args...)

ifft(      A::AffineOperator, args...) = compose(IDFT,       A, args...)

dct(       A::AffineOperator, args...) = compose(DCT,        A, args...)

idct(      A::AffineOperator, args...) = compose(IDCT,       A, args...)

finitediff(A::AffineOperator, args...) = compose(FiniteDiff, A, args...)

variation( A::AffineOperator, args...) = compose(Variation,  A, args...)

conv(      A::AffineOperator, args...) = compose(Conv,       A, args...)

xcorr(     A::AffineOperator, args...) = compose(Xcorr,      A, args...)
