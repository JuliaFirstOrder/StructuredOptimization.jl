immutable Scale{T <: Number} <: LinearOperator
  coeff::T
  A::LinearOperator
end

# So we avoid nesting stuff
Scale{T <: Number}(coeff::T, L::Scale) = Scale(coeff.*L.coeff, L.A);

size(L::Scale) = size(L.A)

  domainType(  L::Scale) =   domainType(L.A)
codomainType(  L::Scale) = codomainType(L.A)
isEye(         L::Scale) = isEye(L.A) 
isDiagonal(    L::Scale) = isDiagonal(L.A) 
isGramDiagonal(L::Scale) = isGramDiagonal(L.A)
isInvertible(  L::Scale) = isInvertible(L.A)

transpose(L::Scale) = Scale(conj(L.coeff),L.A')

function A_mul_B!(y, L::Scale, x)
  A_mul_B!(y, L.A, x)
  y .*= L.coeff
end

# redefine scalar multiplication for convenience
(*){T <: Real   }(c::T, L::LinearOperator) = Scale(c, L)
(*){T <: Complex}(c::T, L::LinearOperator) = 
codomainType(L) <: Complex && domainType(L) <: Complex ? Scale(c, L) : 
error("cannot scale Real operator with Complex scalar")

# redefine unary `-` for convenience
(-)(L::LinearOperator) = Scale(-one(Real), L)

fun_name(L::Scale)  = " $(round(L.coeff,4)) * $(fun_name(L.A))"
fun_type(L::Scale)  = fun_type(L.A)
