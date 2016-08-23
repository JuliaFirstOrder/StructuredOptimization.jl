# indicator of the ball of matrices with (at most) a given rank

"""
  indBallRank(r::Int64)

Returns the function `g = ind{X : rank(X) ⩽ r}`, for an integer parameter `r > 0`.
"""

immutable indBallRank <: ProximableFunction
  r::Int64
  indBallRank(r::Int64) =
    r <= 0 ? error("parameter r must be a positive integer") : new(r)
end

function call(f::indBallRank, x::RealOrComplexMatrix)
  maxr = minimum(size(x))
  if maxr <= f.r return 0.0 end
  u, s, v = svds(x, nsv=f.r+1)
  # the tolerance in the following line should be customizable
  if s[end]/s[1] > 1e-15 return +Inf end
  return 0.0
end

function prox(f::indBallRank, gamma::Float64, x::RealOrComplexMatrix)
  maxr = minimum(size(x))
  if maxr <= f.r return (x, 0.0) end
  u, s, v = svds(x, nsv=f.r)
  return (u*diagm(s))*v', 0.0
end

fun_name(f::indBallRank) = "indicator of the set of rank-r matrices"
fun_type(f::indBallRank) = "C^{n×m} → R ∪ {+∞}"
fun_expr(f::indBallRank) = "x ↦ 0 if rank(x) ⩽ r, +∞ otherwise"
fun_params(f::indBallRank) = "r = $(f.r)"
