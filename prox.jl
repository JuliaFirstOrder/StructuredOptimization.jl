# ------------------------------------------------------------------------------
# prox.jl - library of nonsmooth functions and associated proximal mappings
# ------------------------------------------------------------------------------

# L2 norm (times a constant, or weighted)

function l2norm(lambda::Float64=1.0)
  function prox_l2norm(x, gamma::Float64=1.0)
    y = x/(1+lambda*gamma)
    return y, (lambda/2)*dot(y[:],y[:])
  end
end

function l2norm(lambda::Array{Float64})
  function prox_l2norm(x, gamma::Float64=1.0)
    y = x./(1+lambda*gamma)
    return y, 0.5*dot(lambda[:].*y[:],y[:])
  end
end

# L1 norm (times a constant, or weighted)

function l1norm(lambda::Float64=1.0)
  function prox_l1norm(x, gamma::Float64=1.0)
    y = sign(x).*max(0.0, abs(x)-gamma*lambda)
    return y, lambda*norm(y,1)
  end
end

function l1norm(lambda::Array{Float64})
  function prox_l1norm(x, gamma::Float64=1.0)
    y = sign(x).*max(0.0, abs(x)-gamma*lambda)
    return y, norm(lambda.*y,1)
  end
end

# L0 pseudo-norm (times a constant)

function l0norm(lambda::Float64=1.0)
  function prox_l0norm(x, gamma::Float64=1.0)
    over = abs(x) > sqrt(2*gamma*lambda);
    y = x.*over;
    return y, lambda*countnz(y)
  end
end

# indicator of the L0 norm ball with given (integer) radius

function l0ball(r::Int64)
  function proj_l0ball(x, gamma::Float64=1.0)
    y = zeros(size(x))
    if r < log2(length(x))
      p = selectperm(abs(x), 1:r, rev=true)
      y[p] = x[p]
    else
      p = sortperm(abs(x), rev=true)
      y[p[1:r]] = x[p[1:r]]
    end
    return y, 0.0
  end
end

# indicator of the ball of matrices with (at most) a given rank

function rankball(r::Int64)
  function proj_rankball(x, gamma::Float64=1.0)
    u, s, v = svds(x, nsv=r)
    return (u*diagm(s))*v', 0.0
  end
end
