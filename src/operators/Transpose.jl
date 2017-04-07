# Should we call it "Adjoint"?
immutable Transpose <: LinearOperator
	A::LinearOperator
end

# So we avoid nesting stuff
Transpose(L::Transpose) = L.A;

size(L::Transpose) = size(L.A,2), size(L.A,1)

  domainType(L::Transpose) = codomainType(L.A)
codomainType(L::Transpose) =   domainType(L.A)

A_mul_B!(y, L::Transpose, x) = Ac_mul_B!(y, L.A, x)
Ac_mul_B!(y, L::Transpose, x) = A_mul_B!(y, L.A, x)

fun_name(L::Transpose)  = "$(fun_name(L.A)) (adjoint)"

# redefine `transpose` for convenience
transpose(L::LinearOperator) = Transpose(L)
