export ls

immutable LinearLeastSquares <: SmoothFunction
	lambda::Real
end

function (f::LinearLeastSquares)(x::AbstractArray)
	return 0.5*f.lambda*vecnorm(x)^2
end

function gradient!(grad::AbstractArray, f::LinearLeastSquares)  
	f.lambda == 1. ? nothing : grad .= f.lambda.*grad
end

function get_prox(T::LinearLeastSquares)
	return ProximalOperators.SqrNormL2(T.lambda)
end

ls(x::Variable)       = ls(eye(x))
ls(A::AffineOperator) = CostFunction(variable(A), [LinearLeastSquares(1.)], [A])

*(lambda::Real,f::LinearLeastSquares) = LinearLeastSquares(f.lambda*lambda)

fun_name(f::LinearLeastSquares,i::Int64) = " λ$i/2 ‖A$(i)x‖² "
fun_par( f::LinearLeastSquares,i::Int64)  = " λ$i = $(round(f.lambda,3)) "


