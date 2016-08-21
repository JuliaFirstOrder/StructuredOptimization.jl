# L2,1 norm/Sum of norms of columns or rows (times a constant)

"""
  normL21(λ::Float64=1.0, dim == 1)

Returns the function `g(X) = λsqrt(sum(X.*X, dim))` for a nonnegative parameter
`λ ⩾ 0`, i.e., the (scaled) sum of the norms of the columns (or rows) of X.

If dim is not specified or dim = 1, then the norms are computed column-wise. If
instead dim = 2 then the norms are computed row-wise.
"""

immutable normL21 <: ProximableFunction
  lambda::Float64
  dim::Int
  normL21(lambda::Float64=1.0, dim=1) = lambda < 0 ? error("parameter λ must be nonnegative") : new(lambda, dim)
end

function call(f::normL21, X::RealOrComplexArray)
  return f.lambda*sum(sqrt(sum(abs(X).^2, f.dim)))
end

function prox(f::normL21, gamma::Float64, X::RealOrComplexArray)
  Y = max(0, 1-f.lambda*gamma./sqrt(sum(abs(X).^2, f.dim))).*X
  return Y, f.lambda*sum(sqrt(sum(abs(Y).^2, f.dim)))
end
