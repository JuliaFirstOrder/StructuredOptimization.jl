# indicator of an affine set

"""
  indAffine(A, b)

Returns the function `g = ind{x : Ax = b}`.
"""

immutable indAffine <: ProximableFunction
  A::Array{Float64,2}
  b::Array{Float64,1}
  L::Array{Float64,2}
  function indAffine(A::Array{Float64,2}, b::Array{Float64,1})
    if size(A,1) > size(A,2)
      error("A must be full row rank")
    end
    AAt = A*A';
    L = chol(AAt, Val{:L})
    new(A, b, L)
  end
end

function call(f::indAffine, x::Array{Float64,1})
  if norm(f.A*x - f.b, Inf) > 1e-14 return Inf end
  return 0.0
end

function prox(f::indAffine, gamma::Float64, x::Array{Float64,1})
  res = f.A*x - f.b
  y = x - f.A'*(f.L'\(f.L\res))
  return y, 0.0
end
