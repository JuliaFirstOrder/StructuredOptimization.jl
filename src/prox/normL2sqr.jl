# squared L2 norm (times a constant, or weighted)

immutable normL2sqr{T <: Union{Float64,Array{Float64}}} <: ProximableFunction
  lambda::T
  normL2sqr(lambda::Float64) = lambda < 0 ? error("parameter λ must be nonnegative") : new(lambda)
  normL2sqr(lambda::Array{Float64}) = any(lambda .< 0) ? error("coefficients in λ must be nonnegative") : new(lambda)
end

"""
  normL2sqr(λ::Float64=1.0)

Returns the function `g(x) = (λ/2)(x'x)`, for a real parameter `λ ⩾ 0`.
"""

normL2sqr(lambda::Float64=1.0) = normL2sqr{Float64}(lambda)

"""
  normL2sqr(λ::Array{Float64})

Returns the function `g(x) = (1/2)(λ.*x)'x`, for an array of real parameters `λ ⩾ 0`.
"""

normL2sqr(lambda::Array{Float64}) = normL2sqr{Array{Float64}}(lambda)

function call(f::normL2sqr{Float64}, x::Array{Float64})
  return (f.lambda/2)*vecdot(x,x)
end

function call(f::normL2sqr, x::Array{Float64})
  return 0.5*vecdot(f.lambda.*x,x)
end

function prox(f::normL2sqr{Float64}, gamma::Float64, x::Array{Float64})
  y = x/(1+f.lambda*gamma)
  return y, (f.lambda/2)*vecdot(y,y)
end

function prox(f::normL2sqr, gamma::Float64, x::Array{Float64})
  y = x./(1+f.lambda*gamma)
  return y, 0.5*vecdot(f.lambda.*y,y)
end

fun_type(f::normL2sqr) = "R^n → R"
fun_expr(f::normL2sqr{Float64}) = "x ↦ (λ/2)||x||^2"
fun_expr(f::normL2sqr{Array{Float64}}) = "x ↦ (1/2)sum( λ_i (x_i)^2 )"
fun_params(f::normL2sqr{Float64}) = "λ = $(f.lambda)"
fun_params(f::normL2sqr{Array{Float64}}) = string("λ = ", typeof(f.lambda), " of size ", size(f.lambda))
