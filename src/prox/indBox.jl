# indicator of a generic box

"""
  indBox(lb, ub)

Returns the function `g = ind{x : lb ⩽ x ⩽ ub}`. Parameters `lb` and `ub` can be
either scalars or arrays of the same dimension as `x`, and must satisfy `lb <= ub`.
Bounds are allowed to take values `-Inf` and `+Inf`.
"""

immutable indBox <: ProximableFunction
  lb
  ub
  indBox(lb,ub) = any(lb .> ub) ? error("arguments lb, ub must satisfy lb <= ub") : new(lb, ub)
end

function call(f::indBox, x::Array{Float64})
  if any(x .< f.lb) || any(x .> f.ub) return +Inf end
  return 0.0
end

function prox(f::indBox, gamma::Float64, x::Array{Float64})
  y = min(f.ub, max(f.lb, x))
  return y, 0.0
end

# indicator of the L-infinity ball (box centered in the origin)

"""
  indBallInf(r::Float64)

Returns the indicator function of an infinity-norm ball, that is function
`g(x) = ind{maximum(abs(x)) ⩽ r}` for `r ⩾ 0`.
"""

indBallInf(r::Float64) = indBox(-r, r)

# indicator of the nonnegative orthant

"""
  indNonnegative()

Returns the indicator function the nonnegative orthant, that is

  `g(x) = 0 if x ⩾ 0, +∞ otherwise`
"""

indNonnegative() = indBox(0, +Inf)

fun_name(f::indBox) = "indicator of a box"
fun_type(f::indBox) = "R^n → R ∪ {+∞}"
fun_expr(f::indBox) = "x ↦ 0 if all(lb ⩽ x ⩽ ub), +∞ otherwise"
fun_params(f::indBox) =
  string( "lb = ", typeof(f.lb) <: Array ? string(typeof(f.lb), " of size ", size(f.lb)) : f.lb, ", ",
          "ub = ", typeof(f.ub) <: Array ? string(typeof(f.ub), " of size ", size(f.ub)) : f.ub)
