export ls

immutable LinearLeastSquares <: QuadraticFunction
	lambda::Real
end

lambda(f::LinearLeastSquares) = f.lambda

function (f::LinearLeastSquares)(x::AbstractArray)
	return 0.5*f.lambda*vecnorm(x)^2
end

function gradient!(grad::AbstractArray, f::LinearLeastSquares, x::AbstractArray)  
	grad .= (*).(f.lambda, x)
	return f(x)
end

function gradstep!(x1::AbstractArray, f::LinearLeastSquares, x0::AbstractArray, gamma::Real)
	fx1 = 0.0
	s = (1-gamma*f.lambda)
	for i in eachidex(x)
		x1[i] = s*x0[i]  
		fx1 += x1[i]^2
	end
	return 0.5*f.lambda*fx1
end

function get_prox(T::LinearLeastSquares)
	return ProximalOperators.SqrNormL2(T.lambda)
end

ls(x::Variable)       = ls(eye(x))
ls(A::AffineOperator) = CostFunction(variable(A), [LinearLeastSquares(1.)], [A])

*(lambda::Real,f::LinearLeastSquares) = LinearLeastSquares(f.lambda*lambda)

fun_name(f::LinearLeastSquares,i::Int64) = " λ$i/2 ‖A$(i)x‖² "
fun_par( f::LinearLeastSquares,i::Int64)  = " λ$i = $(round(f.lambda,3)) "


