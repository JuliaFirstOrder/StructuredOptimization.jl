# indicator function of the matrices with at most

"""
  indBallL20(r::Int64, dim=1)

For an integer parameter `r > 0`, if dim=1 then returns the function
`g = ind{X : countnz(||X(:,i)||_2) ⩽ r}`, if dim=2 then instead
`g = ind{X : countnz(||X(i,:)||_2) ⩽ r}`.
"""

immutable indBallL20
  r::Int64
  dim::Int
end

indBallL20(r::Int64) = indBallL20(r, 1)

function call(f::indBallL20, X::RealOrComplexArray)
  if countnz(sqrt(sum(abs(x).^2,dim))) > r return +Inf end
  return 0.0
end

function prox(f::indBallL20, gamma::Float64, X::RealOrComplexArray)
  y = zeros(x)
  if r < log2(size(x,dim))
    p = selectperm( sqrt(sum(abs(x).^2,dim)[:]) , 1:r, rev=true)
    if dim == 1
    y[:,p] = x[:,p]
    elseif dim == 2
    y[p,:] = x[p,:]
    end
  else
    p = sortperm(sqrt(sum(abs(x).^2,dim)[:]), rev=true)
    if dim == 1
    y[:,p[1:r]] = x[:,p[1:r]]
    elseif dim == 2
    y[p[1:r],:] = x[p[1:r],:]
    end
  end
  return y, 0.0
end
