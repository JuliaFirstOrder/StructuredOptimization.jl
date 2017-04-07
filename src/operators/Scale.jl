immutable Scale{T <: Real} <: LinearOperator
  coeff::T
	A::LinearOperator
end

# So we avoid nesting stuff
Scale{T <: Real}(coeff::T, L::Scale) = Scale(coeff*L.coeff, L.A);

size(L::Scale) = size(L.A)

function A_mul_B!(y, L::Scale, x)
  A_mul_B!(y, L.A, x)
  y .*= L.coeff
end

function At_mul_B!(y, A::Scale, x)
  At_mul_B!(y, L.A, x)
  y .*= L.coeff
end

fun_name(L::Scale)  = "$(fun_name(L.A)) (scaled by $(L.coeff))"

# redefine scalar multiplication for convenience
(*){T <: Real}(c::T, L::LinearOperator) = Scale(c, L)

# redefine unary `-` for convenience
(-)(L::LinearOperator) = Scale(-one(Real), L)
