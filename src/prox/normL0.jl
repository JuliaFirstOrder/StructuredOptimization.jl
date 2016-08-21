# L0 pseudo-norm (times a constant)

"""
  normL0(λ::Float64=1.0)

Returns the function `g(x) = λ*countnz(x)`, for a nonnegative parameter `λ ⩾ 0`.
"""

immutable normL0 <: ProximableFunction
  lambda::Float64
  normL0(lambda::Float64=1.0) = lambda < 0 ? error("parameter λ must be nonnegative") : new(lambda)
end

function call(f::normL0, x::Array{Float64})
  return f.lambda*countnz(x)
end

function prox(f::normL0, gamma::Float64, x::Array{Float64})
  over = abs(x) .> sqrt(2*gamma*f.lambda);
  y = x.*over;
  return y, f.lambda*countnz(y)
end
