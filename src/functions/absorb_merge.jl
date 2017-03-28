
# absorb linear operator into proximable operator
absorbOp(A::AffineOperator, p::ProximableFunction) = absorbOp(operator(A), p, tilt(A))

absorbOp{L <:IdentityOperator}(A::L, p::ProximableFunction, b) = b == 0. ? p : Precompose(p, 1., b)
absorbOp{L <:DiagonalOperator}(A::L, p::ProximableFunction, b) = Precompose(p, A.d, b)

# merge Proximal operators 
mergeProx{T<:AffineOperator}(p::ProximableFunction, lambda, A::T) =  Regularize(p, lambda, -tilt(A.b))
mergeProx{T<:LinearOperator}(p::ProximableFunction, lambda, A::T) =  Regularize(p, lambda,     0.    )





