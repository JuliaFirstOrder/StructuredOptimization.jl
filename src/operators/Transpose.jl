# Should we call it "Adjoint"?
immutable Transpose <: LinearOperator
	A::LinearOperator
end

# So we avoid nesting stuff
Transpose(L::Transpose) = L.A;

size(L::Transpose) = size(L.A,2), size(L.A,1)

A_mul_B!(y, L::Transpose, x) = At_mul_B!(y, L.A, x)
At_mul_B!(y, L::Transpose, x) = A_mul_B!(y, L.A, x)

fun_name(L::Transpose)  = "$(fun_name(L.A)) (transposed)"

# redefine `transpose` for convenience
transpose(L::LinearOperator) = Transpose(L)
