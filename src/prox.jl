# ------------------------------------------------------------------------------
# prox.jl - library of nonsmooth functions and associated proximal mappings
# ------------------------------------------------------------------------------

# L2 norm (times a constant, or weighted)

function normL2(lambda::Float64=1.0)
  function call_normL2(x::Array{Float64})
    return (lambda/2)*vecdot(x,x)
  end
  function call_normL2(x::Array{Float64}, gamma::Float64)
    y = x/(1+lambda*gamma)
    return y, (lambda/2)*vecdot(y,y)
  end
  return call_normL2
end

function normL2(lambda::Array{Float64})
  function call_normL2(x::Array{Float64})
    return 0.5*vecdot(lambda.*x,x)
  end
  function call_normL2(x::Array{Float64}, gamma::Float64=1.0)
    y = x./(1+lambda*gamma)
    return y, 0.5*vecdot(lambda.*y,y)
  end
  return call_normL2
end

# L1 norm (times a constant, or weighted)

"""
  normL1(λ::Float64=1.0)

Returns the function `g(x) = λ||x||_1`, for a real parameter `λ ⩾ 0`.
"""
function normL1(lambda::Float64=1.0)
  if lambda < 0 error("parameter λ must be nonnegative") end
  function call_normL1(x::Array{Float64})
    return lambda*vecnorm(x,1)
  end
  function call_normL1(x::Array{Float64}, gamma::Float64)
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

function normL21(lambda::Float64=1.0, dim = 1)
	function prox_l21group(x, gamma::Float64=1.0)
		n = size(x,dim)
		y = max(0, 1-lambda*gamma./sqrt(sum(abs(x).^2,dim))).*x
		return y, lambda*norm(sqrt(sum(abs(y).^2,dim)'),1)
	end
end

# L0 pseudo-norm (times a constant)

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

Returns the function `ind{x : countnz(x) ⩽ r}`, for an integer parameter `r > 0`.
"""
function indBallL0(r::Int64)
  if r < 0 error("parameter r must be positive") end
  function call_indBallL0(x::Array{Float64})
    if countnz(x) > r return +Inf end
    return 0.0
  end
  function call_indBallL0(x::Array{Float64}, gamma::Float64)
    y = zeros(size(x))
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

function indBallRank(r::Int64)
  function call_indBallRank(x::Array{Float64})
    u, s, v = svds(x, nsv=r+1)
    if s[end] > 0 return +Inf end
    return 0.0
  end
  function call_indBallRank(x::Array{Float64}, gamma::Float64)
    u, s, v = svds(x, nsv=r)
    return (u*diagm(s))*v', 0.0
  end
  return call_indBallRank
end

# indicator of a generic box

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

function indBallInf(r::Float64)
  indBox(-r, r)
end
