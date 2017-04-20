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
       
(*)(L::LinearOperator, x::Variable) = LinearTerm(x,L)

(*)(L::LinearOperator,A::LinearTerm)         = LinearTerm(variable(A),L*operator(A))
(*)(L::LinearOperator,A::AffineTerm)   = LinearTerm(variable(A),L*operator(A))+L*tilt(A)

(*){T<:Number}(coeff::T,A::LinearTerm)       = LinearTerm(variable(A),coeff*operator(A))
(*){T<:Number}(coeff::T,A::AffineTerm) = LinearTerm(variable(A),coeff*operator(A))+coeff.*tilt(A)

(*){T<:Number}(coeff::T,x::Variable)     = LinearTerm(x,coeff*Eye(size(x)))
	     
zeros(x::AbstractVariable) = Zeros(~x)*x

eye(x::AbstractVariable) = Eye(~x)*x

reshape(x::AbstractVariable, args...) = Reshape(Eye(~x), args...)*x
reshape(A::LinearTerm          , args...) = LinearTerm(variable(A), Reshape(operator(A), args...))
reshape(A::AffineTerm    , args...) = 
LinearTerm(variable(A), Reshape(operator(A), args...))+reshape(tilt(A), args...)

(.*)(d::AbstractVector,x::AbstractVariable) =   DiagOp(~x, d)*x

getindex(x::AbstractVariable, args...) = GetIndex(~x, args)*x

(*)(M::AbstractMatrix,x::AbstractVariable)  = MatrixOp(~x, M)*x

fft(x::AbstractVariable) = DFT(~x)*x

ifft(x::AbstractVariable) = IDFT(~x)*x

dct(x::AbstractVariable) = DCT(~x)*x

idct(x::AbstractVariable) = IDCT(~x)*x

finitediff(x::AbstractVariable, args...) = FiniteDiff(~x, args...)*x

variation(x::AbstractVariable, args...) = Variation(~x, args...)*x

conv(x::AbstractVariable, args...)  = Conv(~x, args...)*x

xcorr(x::AbstractVariable, args...) = Xcorr(~x, args...)*x


