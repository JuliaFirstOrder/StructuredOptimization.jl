export LinearExpression

immutable LinearExpression{T <: LinearOperator} <: AbstractAffineExpression
	L::T
	x::Variable
	function LinearExpression{T}(L::T, x::Variable) where {T <: LinearOperator}
		if size(L, 2) != size(x)
			throw(ArgumentError("Size of the operator domain $(size(L, 2)) must match size of the variable $(size(x))"))
		end
		if domainType(L) != eltype(x)
			throw(ArgumentError("Type of the operator domain $(domainType(L)) must match type of the variable $(eltype(x))"))
		end
		new(L, x)
	end
end

LinearExpression{T <: LinearOperator}(L::T, x::Variable) = LinearExpression{T}(L, x)

LinearExpression{T <: LinearOperator, E <: LinearExpression}(L::T, e::E) = LinearExpression(Compose(L, e.L), e.x)

# Properties

is_eye(L::LinearExpression) = is_eye(L.L)
is_null(L::LinearExpression) = is_null(L.L)
is_diagonal(L::LinearExpression) = is_diagonal(L.L)
is_gram_diagonal(L::LinearExpression) = is_gram_diagonal(L.L)
is_invertible(L::LinearExpression) = is_invertible(L.L)
is_full_row_rank(L::LinearExpression) = is_full_row_rank(L.L)
is_full_column_rank(L::LinearExpression) = is_full_column_rank(L.L)

# Syntax

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
