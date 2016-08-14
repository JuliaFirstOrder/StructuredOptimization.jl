# ------------------------------------------------------------------------------
# prox.jl - library of nonsmooth functions and associated proximal mappings
# ------------------------------------------------------------------------------

RealOrComplexArray = Union{Array{Float64}, Array{Complex{Float64}}}

# L2 norm (times a constant)

"""
  normL2(λ::Float64=1.0)

Returns the function `g(x) = λ||x||_2`, for a real parameter `λ ⩾ 0`.
"""
function normL2(lambda::Float64=1.0)
  function call_normL2(x::Array{Float64})
    return lambda*vecnorm(x)
  end
  function call_normL2(x::Array{Float64}, gamma::Float64)
    vecnormx = vecnorm(x)
    scale = max(0, 1-lambda*gamma/vecnormx)
    y = scale*x
    return y, lambda*scale*vecnormx
  end
  return call_normL2
end

# squared L2 norm (times a constant, or weighted)

"""
  normL2sqr(λ::Float64=1.0)

Returns the function `g(x) = (λ/2)(x'x)`, for a real parameter `λ ⩾ 0`.
"""
function normL2sqr(lambda::Array{Float64})
  function call_normL2sqr(x::Array{Float64})
    return (lambda/2)*vecdot(x,x)
  end
  function call_normL2sqr(x::Array{Float64}, gamma::Float64=1.0)
    y = x/(1+lambda*gamma)
    return y, (lambda/2)*vecdot(y,y)
  end
  return call_normL2sqr
end

"""
  normL2sqr(λ::Array{Float64})

Returns the function `g(x) = (1/2)(λ.*x)'x`, for an array of real parameters `λ ⩾ 0`.
"""
function normL2sqr(lambda::Array{Float64})
  function call_normL2sqr(x::Array{Float64})
    return 0.5*vecdot(lambda.*x,x)
  end
  function call_normL2sqr(x::Array{Float64}, gamma::Float64=1.0)
    y = x./(1+lambda*gamma)
    return y, 0.5*vecdot(lambda.*y,y)
  end
  return call_normL2sqr
end

# L1 norm (times a constant, or weighted)

"""
  normL1(λ::Float64=1.0)

Returns the function `g(x) = λ||x||_1`, for a real parameter `λ ⩾ 0`.
"""
function normL1(lambda::Float64=1.0)
  if lambda < 0 error("parameter λ must be nonnegative") end
  function call_normL1(x::RealOrComplexArray)
    return lambda*vecnorm(x,1)
  end
  function call_normL1(x::RealOrComplexArray, gamma::Float64)
	  y = sign(x).*max(0.0, abs(x)-gamma*lambda)
	  return y, lambda*vecnorm(y,1)
  end
  return call_normL1
end

"""
  normL1(λ::Array{Float64})

Returns the function `g(x) = sum(λ_i|x_i|, i = 1,...,n)`, for a vector of real
parameters `λ_i ⩾ 0`.
"""
function normL1(lambda::Array{Float64})
  if any(lambda .< 0) error("coefficients in λ must be nonnegative") end
  function call_normL1(x::Array{Float64})
    return vecnorm(lambda.*x,1)
  end
  function call_normL1(x::Array{Float64}, gamma::Float64)
    y = sign(x).*max(0.0, abs(x)-gamma*lambda)
    return y, vecnorm(lambda.*y,1)
  end
  return call_normL1
end

# L2,1 norm/Sum of norms (times a constant)

function normL21(lambda::Float64=1.0, dim=1)
	function call_normL21(x::Array{Float64,2}, gamma::Float64)
		y = max(0, 1-lambda*gamma./sqrt(sum(x.^2, dim))).*x
		return y, lambda*norm(sqrt(sum(y.^2, dim)'),1)
	end
end

# L0 pseudo-norm (times a constant)

"""
  normL0(λ::Float64=1.0)

Returns the function `g(x) = λ*countnz(x)`, for a nonnegative parameter `λ ⩾ 0`.
"""
function normL0(lambda::Float64=1.0)
  function call_normL0(x::Array{Float64})
    return lambda*countnz(x)
  end
  function call_normL0(x::Array{Float64}, gamma::Float64)
    over = abs(x) .> sqrt(2*gamma*lambda);
    y = x.*over;
    return y, lambda*countnz(y)
  end
  return call_normL0
end

# indicator of the L0 norm ball with given (integer) radius

"""
  indBallL0(r::Int64)

Returns the function `g(x) = ind{countnz(x) ⩽ r}`, for an integer parameter `r > 0`.
"""
function indBallL0(r::Int64)
  if r < 0 error("parameter r must be positive") end
  function call_indBallL0(x::RealOrComplexArray)
    if countnz(x) > r return +Inf end
    return 0.0
  end
  function call_indBallL0(x::RealOrComplexArray, gamma::Float64)
    y = zeros(x)
    if r < log2(length(x))
      p = selectperm(abs(x)[:], 1:r, rev=true)
      y[p] = x[p]
    else
      p = sortperm(abs(x)[:], rev=true)
      y[p[1:r]] = x[p[1:r]]
    end
    return y, 0.0
  end
  return call_indBallL0
end

# indicator of the ball of matrices with (at most) a given rank

"""
  indBallRank(r::Int64)

Returns the function `g(X) = ind{rank(X) ⩽ r}`, for an integer parameter `r > 0`.
"""
function indBallRank(r::Int64)
  function call_indBallRank(x::RealOrComplexArray)
    u, s, v = svds(x, nsv=r+1)
    if s[end]/s[1] >= 1e-15 return +Inf end
    return 0.0
  end
  function call_indBallRank(x::RealOrComplexArray, gamma::Float64)
    u, s, v = svds(x, nsv=r)
    return (u*diagm(s))*v', 0.0
  end
  return call_indBallRank
end

# indicator of a generic box

"""
  indBox(lb, ub)

Returns the function `g(x) = ind{lb ⩽ x ⩽ ub}`. Parameters `lb` and `ub` can be
either scalars or arrays of the same dimension as `x`.
"""
function indBox(lb, ub)
  function call_indBox(x::Array{Float64})
    if any(x .< lb) || any(x .> ub) return +Inf end
    return 0.0
  end
  function call_indBox(x::Array{Float64}, gamma::Float64)
    y = min(ub, max(lb, x))
    return y, 0.0
  end
  return call_indBox
end

# indicator of the L-infinity ball (box centered in the origin)

"""
  indBallInf(r::Float64)

Returns the indicator function of an infinity-norm ball, that is function
`g(x) = ind{maximum(abs(x)) ⩽ r}` for `r ⩾ 0`.
"""
function indBallInf(r::Float64)
  indBox(-r, r)
end
