export ls

type LinearLeastSquares{T <: AffineOperator} <: QuadraticTerm
	A::T
	At::LinearOperator
	lambda::Number
end

function (t::LinearLeastSquares)(x::AbstractArray)
	0.5*t.lambda*deepvecnorm(x)^2
end

ls{T <: AffineOperator}(A::T) = LinearLeastSquares(A, A',1)
ls(x::OptVar) = LinearLeastSquares(eye(x),eye(x),1)

function get_prox(T::LinearLeastSquares)
	return ProximalOperators.SqrNormL2(T.lambda)
end

*(lambda::Number, T::LinearLeastSquares) = LinearLeastSquares(T.A,lambda)

fun_name(T::LinearLeastSquares) = " λ/2 ‖⋅‖² "
fun_par(T::LinearLeastSquares)  = " λ = $(round(T.lambda,3)) "


function gradient!(grad::AbstractArray, t::LinearLeastSquares, x::AbstractArray)
	A_mul_B!(grad,t.At,x)
	grad .*= t.lambda
end

gradient(t::LinearLeastSquares, x::AbstractArray) = t.lambda.*(t.A'*x)

function evaluate!(resx::AbstractArray, t::LinearLeastSquares, x::AbstractArray)
	A_mul_B!(resx,t.A,x)
	return t(resx)
end

function evaluate(t::LinearLeastSquares, x::AbstractArray)
	resx = t.A*x
	return resx, t(resx)
end
