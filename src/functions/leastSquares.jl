immutable LeastSquares <: QuadraticFunction
	lambda::Real
end

lambda(f::LeastSquares) = f.lambda

function (f::LeastSquares)(x::AbstractArray)
	return 0.5*f.lambda*vecnorm(x)^2
end

function gradient!(grad::AbstractArray, f::LeastSquares, x::AbstractArray)
	grad .= (*).(f.lambda, x)
	return f(x)
end

function gradstep!(x1::AbstractArray, f::LeastSquares, x0::AbstractArray, gamma::Real)
	fx1 = 0.0
	s = (1-gamma*f.lambda)
	for i in eachindex(x0)
		x1[i] = s*x0[i]
		fx1 += x1[i]^2
	end
	return 0.5*f.lambda*fx1
end

function get_prox(T::LeastSquares)
	return ProximalOperators.SqrNormL2(T.lambda)
end

*(lambda::Real, f::LeastSquares) = LeastSquares(f.lambda*lambda)

fun_name(f::LeastSquares,i::Int64) = " λ$i/2 ‖A$(i)x‖² "
fun_par( f::LeastSquares,i::Int64)  = " λ$i = $(round(f.lambda,3)) "
