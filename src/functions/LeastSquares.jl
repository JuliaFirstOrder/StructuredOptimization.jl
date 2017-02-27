type LeastSquares{T <: AffineOperator} <: SmoothTerm
	A::T
	lambda::Number
end

ls{T <: AffineOperator}(A::T) = LeastSquares(A,1)
ls(x::OptVar) = LeastSquares(eye(x),1)

function get_prox(T::LeastSquares)
	return ProximalOperators.SqrNormL2(T.lambda)
end

*(lambda::Number, T::LeastSquares) = LeastSquares(T.A,lambda)

export ls
