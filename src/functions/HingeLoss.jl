immutable HingeLoss{R1 <: Real, R2 <: Real}   <: NonSmoothFunction
	b::Array{R1,1}
	mu::R2
end

lambda(f::HingeLoss) = f.mu

*(mu::Real,T::HingeLoss) = HingeLoss(T.b,mu)

get_prox(T::HingeLoss) = ProximalOperators.HingeLoss(T.b, T.mu)

fun_name(T::HingeLoss,i::Int64) = " μ$(i) ∑_i max(0, 1 - b$(i)_i A$(i)x )  "
fun_par( T::HingeLoss,i::Int64) = " μ$(i) = $(round(T.mu,3)), b$(i) = $(typeof(T.b)) "
