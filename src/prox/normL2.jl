# L2 norm (times a constant)

"""
  normL2(λ::Float64=1.0)

Returns the function `g(x) = λ||x||_2`, for a real parameter `λ ⩾ 0`.
"""

immutable normL2 <: ProximableFunction
  lambda::Float64
  normL2(lambda::Float64=1.0) = lambda < 0 ? error("parameter λ must be nonnegative") : new(lambda)
end

function call(f::normL2, x::Array{Float64})
  return f.lambda*vecnorm(x)
end

function prox(f::normL2, gamma::Float64, x::Array{Float64})
  vecnormx = vecnorm(x)
  scale = max(0, 1-f.lambda*gamma/vecnormx)
  y = scale*x
  return y, f.lambda*scale*vecnormx
end

fun_name(f::normL2) = "weighted Euclidean norm"
fun_type(f::normL2) = "R^n → R"
fun_expr(f::normL2) = "x ↦ λ||x||_2"
fun_params(f::normL2) = "λ = $(f.lambda)"
