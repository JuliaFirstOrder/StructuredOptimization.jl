# function compose(f::Type, B::LinearExpression, args...)
# 	mid = Array{codomainType(operator(B))}(size(operator(B),1))
# 	A = f(mid, args...)
# 	C = Compose(A,operator(B),mid)
# 	LinearExpression(variable(B),C)
# end
# 
# function compose(f::Type, B::AffineExpression, args...)
# 	mid = Array{codomainType(operator(B))}(size(operator(B),1))
# 	A = f(mid, args...)
# 	C = Compose(A,operator(B),mid)
# 	LinearExpression(variable(B),C)+A*tilt(B)
# end
#
#
# (.*)(d::AbstractVector,A::AbstractAffineExpression  ) =   compose(DiagOp, A, d)
#
# (*)(M::AbstractMatrix,A::AbstractAffineExpression  )  = compose(MatrixOp, A, M)
#
# eye(       A::AbstractAffineExpression         ) = A
#
# zeros(     A::AbstractAffineExpression, args...) = compose(Zeros,      A, args...)
#
# getindex(  A::AbstractAffineExpression, args...) = compose(GetIndex,   A, args)
#
# fft(       A::AbstractAffineExpression, args...) = compose(DFT,        A, args...)
#
# ifft(      A::AbstractAffineExpression, args...) = compose(IDFT,       A, args...)
#
# dct(       A::AbstractAffineExpression, args...) = compose(DCT,        A, args...)
#
# idct(      A::AbstractAffineExpression, args...) = compose(IDCT,       A, args...)
#
# finitediff(A::AbstractAffineExpression, args...) = compose(FiniteDiff, A, args...)
#
# variation( A::AbstractAffineExpression, args...) = compose(Variation,  A, args...)
#
# conv(      A::AbstractAffineExpression, args...) = compose(Conv,       A, args...)
#
# xcorr(     A::AbstractAffineExpression, args...) = compose(Xcorr,      A, args...)
#
# filt(     A::AbstractAffineExpression, args...) = compose(Filt,      A, args...)
#
# mimofilt( A::AbstractAffineExpression, args...) = compose(MIMOFilt,      A, args...)
#
# zeropad(  A::AbstractAffineExpression, args...) = compose(ZeroPad,      A, args...)
