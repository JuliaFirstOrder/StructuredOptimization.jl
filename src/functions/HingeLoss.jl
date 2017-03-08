export hingeloss

type HingeLoss{T <: AffineOperator, R1 <: Real, R2 <: Real}   <: NonSmoothTerm
	A::T
	b::Array{R1,1}
	mu::R2
end

hingeloss(x::OptVar, args...) = hingeloss(eye(x), args...)

hingeloss{T <: AffineOperator, R <: Real}(A::T, b::Array{R,1}) = HingeLoss(A, b, 1.0)
*(mu::Real,T::HingeLoss) = HingeLoss(T.A,T.b,mu)

get_prox(T::HingeLoss) = ProximalOperators.HingeLoss(T.b, T.mu)
