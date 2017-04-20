abstract NormFunction <: NonSmoothFunction

type NormL0  <: NormFunction
	lambda::Real
end

lambda(f::NormFunction) = f.lambda

fun_name(T::NormL0, i::Int64) = " λ$(i) ‖A$(i)x‖₀ "
fun_par(T::NormL0, i::Int64) = " λ$(i) = $(round(T.lambda,3)) "

type IndBallL0 <: NonSmoothFunction
	r::Integer
end

fun_name(T::IndBallL0, i::Int64) = "Ind{‖A$(i)x‖₀ ≤ r$(i)}(x) "
fun_par( T::IndBallL0, i::Int64) = " r$(i) = $(T.r) "

type NormL1  <: NormFunction
	lambda::Real
end

fun_name(T::NormL1, i::Int64) = " λ$(i) ‖A$(i)x‖₁ "
fun_par( T::NormL1, i::Int64) = " λ$(i) = $(round(T.lambda,3)) "

type IndBallL1 <: NonSmoothFunction
	r::Real
end

fun_name(T::IndBallL1, i::Int64) = " Ind{‖A$(i)x‖₁ ≤ r$(i)}(x) "
fun_par( T::IndBallL1, i::Int64) = " r$(i) = $(round(T.r,3)) "

type NormL2 <: NormFunction
	lambda::Real
end

fun_name(T::NormL2, i::Int64) = " λ$(i) ‖A$(i)x‖₂ "
fun_par( T::NormL2, i::Int64) = " λ$(i) = $(round(T.lambda,3)) "

type NormL21 <: NormFunction
	lambda::Real
	dim::Int
end

fun_name(T::NormL21, i::Int64) =  T.dim == 1 ? " λ$(i) ∑_j ‖ (A$(i)x)[:,j] ‖₂ " : 
                                               " λ$(i) ∑_j ‖ (A$(i)x)[j,:] ‖₂ "
fun_par( T::NormL21, i::Int64) = " λ$(i) = $(round(T.lambda,3)) "

type IndBallL2 <: NonSmoothFunction
	r::Real
end

fun_name(T::IndBallL2, i::Int64) = " Ind{‖A$(i)x‖₂ ≤ r$(i)}(x) "
fun_par( T::IndBallL2, i::Int64) = " r$(i) = $(round(T.r,3)) "

type IndSphereL2 <: NonSmoothFunction
	r::Real
end

fun_name(T::IndSphereL2, i::Int64) = " Ind{‖A$(i)x‖₂ = r$(i)}(x) "
fun_par( T::IndSphereL2, i::Int64) = " r$(i) = $(round(T.r,3)) "

type NormLinf <: NormFunction
	lambda::Real
end

fun_name(T::NormLinf, i::Int64) = " λ$(i) ‖A$(i)x‖∞ "
fun_par( T::NormLinf, i::Int64) = " λ$(i) = $(round(T.lambda,3)) "

type IndBallLinf <: NonSmoothFunction
	r::Real
end

fun_name(T::IndBallLinf, i::Int64) = " Ind{‖A$(i)x‖∞ ≤ r$(i)}(x) "
fun_par( T::IndBallLinf, i::Int64) = " r$(i) = $(round(T.r,3)) "

*(lambda::Real,  T::NormL0)   = NormL0(  T.lambda*lambda)
*(lambda::Real,  T::NormL1)   = NormL1(  T.lambda*lambda)
*(lambda::Real,  T::NormL2)   = NormL2(  T.lambda*lambda)
*(lambda::Real,  T::NormL21)  = NormL21( T.lambda*lambda, T.dim)
*(lambda::Real,  T::NormLinf) = NormLinf(T.lambda*lambda)

get_prox(T::NormL0)   = ProximalOperators.NormL0(T.lambda)
get_prox(T::NormL1)   = ProximalOperators.NormL1(T.lambda)
get_prox(T::NormL2)   = ProximalOperators.NormL2(T.lambda)
get_prox(T::NormL21)  = ProximalOperators.NormL21(T.lambda,T.dim)
get_prox(T::NormLinf) = ProximalOperators.NormLinf(T.lambda)

get_prox(T::IndBallL0)   = ProximalOperators.IndBallL0(T.r)
get_prox(T::IndBallL1)   = ProximalOperators.IndBallL1(T.r)
get_prox(T::IndBallL2)   = ProximalOperators.IndBallL2(T.r)
get_prox(T::IndBallLinf) = ProximalOperators.IndBallLinf(T.r)
get_prox(T::IndSphereL2) = ProximalOperators.IndSphereL2(T.r)

