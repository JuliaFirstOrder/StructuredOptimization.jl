# elastic-net regularization

"""
  elasticNet(λ::Float64=1.0, μ::Float64=1.0)

Returns the function `g(x) = μ||x||_1 + (λ/2)||x||^2`, for a real parameters `μ, λ ⩾ 0`.
"""

immutable elasticNet <: ProximableFunction
  mu::Float64
  lambda::Float64
  elasticNet(mu::Float64=1.0, lambda::Float64=1.0) = lambda < 0 || mu < 0 ? error("parameters μ, λ must be nonnegative") : new(mu, lambda)
end

function call(f::elasticNet, x::RealOrComplexArray)
  return f.mu*vecnorm(x,1) + f.lambda*vecnorm(x,2)^2
end

function prox(f::elasticNet, gamma::Float64, x::RealOrComplexArray)
  uz = max(0, abs(x) - gamma*f.mu)/(1 + f.lambda*gamma);
  prox = sign(x).*uz;
  g = f.mu*sum(uz) + (f.lambda/2)*vecnorm(uz)^2;
end
