# These (*) could be all implemented at once using the Postcompose from
# ProximalOperators. Then special cases of Postcompose, where the coefficient
# gets included in f rather than embedding f in an outer object, we should
# probably implement directly in ProximalOperators.

*(lambda::Real, f::NormL0) = NormL0(f.lambda*lambda)
*(lambda::Real, f::NormL1) = NormL1(f.lambda*lambda)
*(lambda::Real, f::NormL2) = NormL2(f.lambda*lambda)
*(lambda::Real, f::NormL21) = NormL21(f.lambda*lambda, f.dim)
*{R <: Real}(lambda::R, f::Conjugate{IndBallL1{R}}) = NormLinf(lambda*f.f.r)
*(lambda::Real, f::Postcompose) = Postcompose(f.f, lambda*f.a, lambda*f.b)

<=(f::NormL0, r::Integer) = IndBallL0(r)
<=(f::NormL1, r::Real) = IndBallL1(r/f.lambda)
<=(f::NormL2, r::Real) = IndBallL2(r/f.lambda)
<={R <: Real}(f::Conjugate{IndBallL1{R}}, r::R) = IndBallLinf(r/f.f.r)

==(f::NormL2, r::Real) = IndSphereL2(r/f.lambda)

sum(f::NormL2, d::Int) = NormL21(f.lambda, d)
