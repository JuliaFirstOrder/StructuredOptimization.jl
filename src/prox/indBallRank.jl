# indicator of the ball of matrices with (at most) a given rank

"""
  indBallRank(r::Int64)

Returns the function `g = ind{X : rank(X) ⩽ r}`, for an integer parameter `r > 0`.
"""

immutable indBallRank <: ProximableFunction
  r::Int64
  indBallRank(r::Int64) = r <= 0 ? error("parameter r must be a positive integer") : new(r)
end

function call(f::indBallRank, x::RealOrComplexArray)
  u, s, v = svds(x, nsv=f.r+1)
  if s[end]/s[1] >= 1e-15 return +Inf end
  return 0.0
end

function prox(f::indBallRank, gamma::Float64, x::RealOrComplexArray)
  u, s, v = svds(x, nsv=f.r)
  return (u*diagm(s))*v', 0.0
end
