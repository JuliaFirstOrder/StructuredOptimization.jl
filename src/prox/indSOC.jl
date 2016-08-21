# indicator of second-order cones

"""
  indSOC(n::Int64)

Returns the indicator of the second-order cone (ice-cream cone) of R^n.
"""

immutable indSOC <: ProximableFunction
  n::Int64
  indSOC(n::Int64) = new(n)
end

function call(f::indSOC, x::Array{Float64,1})
  # the tolerance in the following line should be customizable
  if norm(x[2:end]) - x[1] > 1e-15 return +Inf end
  return 0.0
end

function prox(f::indSOC, gamma::Float64, x::Array{Float64,1})
  nx = norm(x[2:end])
  t = x[1]
  if nx <= -t
    y = zeros(x)
  elseif nx <= t
    y = x
  else
    y = zeros(x)
    r = 0.5 * (1 + t / nx)
    y[1] = r * nx
    y[2:end] = r * x[2:end]
  end
  return y, 0.0
end

fun_name(f::indSOC) = "indicator of a second-order cone"
fun_type(f::indSOC) = "R^n → R ∪ {+∞}"
fun_expr(f::indSOC) = "x ↦ 0 if x[1] >= ||x[2:end]||, +∞ otherwise"
fun_params(f::indSOC) = "n = $(f.n)"
