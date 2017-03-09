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

fun_name(T::HingeLoss) = " μ ∑_i max(0, 1 - b_i ⋅ )  "
fun_par(T::HingeLoss) = " μ = $(round(T.mu,3)), b = $(typeof(T.b)) "
