function compose(f::Type, B::LinearTerm, args...)
	mid = Array{codomainType(operator(B))}(size(operator(B),1))
	A = f(mid, args...)
	C = Compose(A,operator(B),mid)
	LinearTerm(variable(B),C)
end

function compose(f::Type, B::AffineTerm, args...)
	mid = Array{codomainType(operator(B))}(size(operator(B),1))
	A = f(mid, args...)
	C = Compose(A,operator(B),mid)
	LinearTerm(variable(B),C)+A*tilt(B)
end


(.*)(d::AbstractVector,A::AbstractAffineTerm  ) =   compose(DiagOp, A, d)

(*)(M::AbstractMatrix,A::AbstractAffineTerm  )  = compose(MatrixOp, A, M)

eye(       A::AbstractAffineTerm         ) = A

zeros(     A::AbstractAffineTerm, args...) = compose(Zeros,      A, args...)

getindex(  A::AbstractAffineTerm, args...) = compose(GetIndex,   A, args)

fft(       A::AbstractAffineTerm, args...) = compose(DFT,        A, args...)

ifft(      A::AbstractAffineTerm, args...) = compose(IDFT,       A, args...)

dct(       A::AbstractAffineTerm, args...) = compose(DCT,        A, args...)

idct(      A::AbstractAffineTerm, args...) = compose(IDCT,       A, args...)

finitediff(A::AbstractAffineTerm, args...) = compose(FiniteDiff, A, args...)

variation( A::AbstractAffineTerm, args...) = compose(Variation,  A, args...)

conv(      A::AbstractAffineTerm, args...) = compose(Conv,       A, args...)

xcorr(     A::AbstractAffineTerm, args...) = compose(Xcorr,      A, args...)

filt(     A::AbstractAffineTerm, args...) = compose(Filt,      A, args...)

mimofilt( A::AbstractAffineTerm, args...) = compose(MIMOFilt,      A, args...)

zeropad(  A::AbstractAffineTerm, args...) = compose(ZeroPad,      A, args...)
