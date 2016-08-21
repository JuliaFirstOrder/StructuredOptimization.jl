# indicator of the L2 norm ball with given radius

"""
  indBallL2(r::Float64)

Returns the function `g = ind{x : ||x|| ⩽ r}`, for a real parameter `r > 0`.
"""

immutable indBallL2 <: ProximableFunction
  r::Float64
  indBallL2(r::Float64) = r <= 0 ? error("parameter r must be positive") : new(r)
end

function call(f::indBallL2, x::RealOrComplexArray)
  if vecnorm(x) > f.r return +Inf end
  return 0.0
end

function prox(f::indBallL2, gamma::Float64, x::RealOrComplexArray)
  y = x*min(1, f.r/vecnorm(x))
  return y, 0.0
end

fun_name(f::indBallL2) = "indicator of an L2 norm ball"
fun_type(f::indBallL2) = "C^n → R ∪ {+∞}"
fun_expr(f::indBallL2) = "x ↦ 0 if ||x|| ⩽ r, +∞ otherwise"
fun_params(f::indBallL2) = "r = $(f.r)"
