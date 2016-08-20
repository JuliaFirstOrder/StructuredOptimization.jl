# indicator of a generic box

"""
  indBox(lb, ub)

Returns the function `g(x) = ind{lb ⩽ x ⩽ ub}`. Parameters `lb` and `ub` can be
either scalars or arrays of the same dimension as `x`, and must satisfy `lb <= ub`.
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
