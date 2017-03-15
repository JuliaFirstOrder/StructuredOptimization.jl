export ls

type LinearLeastSquares{T <: AffineOperator} <: QuadraticTerm
	A::T
	lambda::Number
end

ls{T <: AffineOperator}(A::T) = LinearLeastSquares(A,1)
ls(x::OptVar) = LinearLeastSquares(eye(x),1)

function get_prox(T::LinearLeastSquares)
	return ProximalOperators.SqrNormL2(T.lambda)
end

*(lambda::Number, T::LinearLeastSquares) = LinearLeastSquares(T.A,lambda)

fun_name(T::LinearLeastSquares) = " λ/2 ‖⋅‖² "
fun_par(T::LinearLeastSquares)  = " λ = $(round(T.lambda,3)) "

