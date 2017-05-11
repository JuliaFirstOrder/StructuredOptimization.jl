immutable Scale{T<:RealOrComplex, 
		C <: AbstractArray{T}, 
		D <:AbstractArray, 
		L <: LinearOperator} <: LinearOperator
  coeff::T
  A::L
end

size(L::Scale) = size(L.A)

# Constructors
Scale{T1 <: Number, T2<:LinearOperator}(coeff::T1, L::T2) = 
Scale{codomainType(L),
      Array{codomainType(L),ndims(L,1)},
      Array{domainType(L),ndims(L,2)},
      typeof(L)}(convert(codomainType(L),coeff), L)
Scale{T <: Number}(coeff::T, L::Scale) = Scale(coeff.*L.coeff, L.A)

## redefine scalar multiplication for convenience
(*){T <: Real   }(c::T, L::LinearOperator) = Scale(convert(codomainType(L),c), L)
(*){T <: Complex}(c::T, L::LinearOperator) = 
codomainType(L) <: Complex ? Scale(convert(codomainType(L),c), L) : 
error("cannot scale Real operator with Complex scalar")

(*){T<:Real   }(coeff::T, L2::DiagOp) = DiagOp(coeff.*L2.d)
(*){T<:Complex}(coeff::T, L2::DiagOp) = DiagOp(coeff.*L2.d)

## redefine unary `-` for convenience
(-)(L::LinearOperator) = Scale(-one(Real), L)

# Operators

function A_mul_B!{T1,C,D,T2<:LinearOperator}(y::C, 
					     L::Scale{T1,C,D,T2}, 
					     x::D)
  A_mul_B!(y, L.A, x)
  y .*= L.coeff
end

# Transformations
transpose(L::Scale) = Scale(conj(L.coeff),L.A')
inv(L::Scale) = Scale(1/L.coeff,inv(L.A))

# Properties
  domainType(  L::Scale) =   domainType(L.A)
codomainType(  L::Scale) = codomainType(L.A)
isEye(         L::Scale) = isEye(L.A) 
isDiagonal(    L::Scale) = isDiagonal(L.A) 
isGramDiagonal(L::Scale) = isGramDiagonal(L.A)
isInvertible(  L::Scale) = isInvertible(L.A)
isScaled(      L::Scale) = true

fun_name(L::Scale)  = " $(round(L.coeff,4)) * $(fun_name(L.A))"
fun_type(L::Scale)  = fun_type(L.A)
