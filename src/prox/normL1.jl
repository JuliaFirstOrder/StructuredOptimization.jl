# L1 norm (times a constant, or weighted)

"""
  normL1(λ::Float64=1.0)

Returns the function `g(x) = λ||x||_1`, for a real parameter `λ ⩾ 0`.

  normL1(λ::Array{Float64})

Returns the function `g(x) = sum(λ_i|x_i|, i = 1,...,n)`, for a vector of real
parameters `λ_i ⩾ 0`.
"""

immutable normL1{T <: Union{Float64,Array{Float64}}} <: ProximableFunction
  lambda::T
end

function call(f::normL1{Float64}, x::RealOrComplexArray)
  return f.lambda*vecnorm(x,1)
end

function call(f::normL1, x::RealOrComplexArray)
  return vecnorm(f.lambda.*x,1)
end

function prox(f::normL1{Float64}, gamma::Float64, x::RealOrComplexArray)
  y = sign(x).*max(0.0, abs(x)-gamma*f.lambda)
  return y, f.lambda*vecnorm(y,1)
end

function prox(f::normL1, gamma::Float64, x::RealOrComplexArray)
  y = sign(x).*max(0.0, abs(x)-gamma*lambda)
  return y, vecnorm(lambda.*y,1)
end
