import Base: reshape,
	     .*,
	     getindex,
	     *,
	     fft,
	     ifft,
	     dct,
	     idct,
	     conv,
	     xcorr,
	     filt

export finitediff,
	     variation,
	     mimofilt,
	     zeropad

# (*)(L::LinearOperator, x) = LinearExpression(L, x)
# (*)(L::LinearOperator, A::LinearExpression) = LinearExpression(variable(A), L*operator(A))
# (*)(L::LinearOperator, A::AffineExpression) = LinearExpression(variable(A), L*operator(A)) + L*tilt(A)
#
# (*){T<:Number}(coeff::T, A::LinearExpression) = LinearExpression(variable(A), coeff*operator(A))
# (*){T<:Number}(coeff::T, A::AffineExpression) = LinearExpression(variable(A), coeff*operator(A)) + coeff.*tilt(A)
#
# (*){T<:Number}(coeff::T, x::Variable) = LinearExpression(x, coeff*Eye(size(x)))

# zeros(x::Variable) = Zeros(~x)*x

# eye(x::Variable) = LinearExpression(Eye(eltype(x), size(x)), x)

# LinearExpression construction

(.*)(d::AbstractVector, x::Variable) =
	LinearExpression(DiagOp(d, eltype(x)), x)

function (*)(M::AbstractMatrix, x::Variable)
	if ndims(x) == 1
		return LinearExpression(MatrixOp(M, eltype(x), 1), x)
	elseif ndims(x) == 2
		return LinearExpression(MatrixOp(M, eltype(x), size(x, 2)), x)
	else
		error("cannot multiply a Matrix by a n-dimensional Variable with n > 2")
	end
end

getindex(x::Variable, args...) = GetIndex(~x, args)*x

fft(x::Variable) = DFT(~x)*x

ifft(x::Variable) = IDFT(~x)*x

dct(x::Variable) = DCT(~x)*x

idct(x::Variable) = IDCT(~x)*x

finitediff(x::Variable, args...) = FiniteDiff(~x, args...)*x

variation(x::Variable, args...) = Variation(~x, args...)*x

conv(x::Variable, args...)  = Conv(~x, args...)*x

xcorr(x::Variable, args...) = Xcorr(~x, args...)*x

filt(x::Variable, args...) = Filt(~x, args...)*x

mimofilt(x::Variable, args...) = MIMOFilt(~x, args...)*x

zeropad(x::Variable, args...) = ZeroPad(~x, args...)*x

reshape(x::Variable, args...) = Reshape(Eye(~x), args...)*x
reshape(A::LinearExpression, args...) = LinearExpression(variable(A), Reshape(operator(A), args...))

# reshape(A::AffineExpression, args...) =
# LinearExpression(variable(A), Reshape(operator(A), args...))+reshape(tilt(A), args...)

# LinearExpression manipulation

import Base: +, -

(+)(ex::LinearExpression) = ex1
(-)(ex::LinearExpression) = LinearExpression(Scale(-1.0, ex.L), ex.x)

# AffineExpression construction

(+)(ex::LinearExpression, b::AbstractArray) = AffineExpression([ex], b)
(-)(ex::LinearExpression, b::AbstractArray) = AffineExpression([ex], -b)
(+)(b::AbstractArray, ex::LinearExpression) = AffineExpression([ex], b)
(-)(b::AbstractArray, ex::LinearExpression) = AffineExpression([-ex], b)

(+)(lex1::LinearExpression, lex2::LinearExpression) = AffineExpression([lex1, lex2])
(-)(lex1::LinearExpression, lex2::LinearExpression) = AffineExpression([lex1, -lex2])

(+)(ex::AffineExpression, b::AbstractArray) = AffineExpression(ex.Ls, ex.b+b)
(-)(ex::AffineExpression, b::AbstractArray) = AffineExpression(ex.Ls, ex.b-b)
(+)(b::AbstractArray, ex::AffineExpression) = AffineExpression(ex.Ls, ex.b+b)
(-)(b::AbstractArray, ex::AffineExpression) = AffineExpression(.-ex.Ls, b-ex.b)

(+)(ex::AffineExpression, lex::LinearExpression) = AffineExpression(vcat(ex.Ls, lex), ex.b)
(-)(ex::AffineExpression, lex::LinearExpression) = AffineExpression(vcat(ex.Ls, -lex), ex.b)
(+)(lex::LinearExpression, ex::AffineExpression) = AffineExpression(vcat(ex.Ls, lex), ex.b)
(-)(lex::LinearExpression, ex::AffineExpression) = AffineExpression(vcat(.-ex.Ls, lex), ex.b)

(+)(ex1::AffineExpression, ex2::AffineExpression) = AffineExpression(vcat(ex1.Ls, ex2.Ls), ex1.b+ex2.b)
(-)(ex1::AffineExpression, ex2::AffineExpression) = AffineExpression(vcat(ex1.Ls, .-ex2.Ls), ex1.b-ex2.b)
