import Base: norm, *, <=, == 

type NormL0{T <: AffineOperator}   <: NonSmoothTerm
	A::T
	lambda::Real
end

fun_name(T::NormL0) = " λ ‖⋅‖₀ "
fun_par(T::NormL0) = " λ = $(round(T.lambda,3)) "

type IndBallL0{T <: AffineOperator} <: NonSmoothTerm
	A::T
	r::Integer
end

fun_name(T::IndBallL0) = " ‖⋅‖₀ ≤ r "
fun_par(T::IndBallL0) = " r = $(T.r) "

type NormL1{T <: AffineOperator}   <: NonSmoothTerm
	A::T
	lambda::Real
end

fun_name(T::NormL1) = " λ ‖⋅‖₁ "
fun_par(T::NormL1) = " λ = $(round(T.lambda,3)) "

type IndBallL1{T <: AffineOperator} <: NonSmoothTerm
	A::T
	r::Real
end

fun_name(T::IndBallL1) = " ‖⋅‖₁ ≤ r "
fun_par(T::IndBallL1) = " r = $(round(T.r,3)) "

type NormL2{T <: AffineOperator} <: NonSmoothTerm
	A::T
	lambda::Real
end

fun_name(T::NormL2) = " λ ‖⋅‖₂ "
fun_par(T::NormL2) = " λ = $(round(T.lambda,3)) "

type IndBallL2{T <: AffineOperator} <: NonSmoothTerm
	A::T
	r::Real
end

fun_name(T::IndBallL2) = " ‖⋅‖₂ ≤ r "
fun_par(T::IndBallL2) = " r = $(round(T.r,3)) "

type IndSphereL2{T <: AffineOperator} <: NonSmoothTerm
	A::T
	r::Real
end

fun_name(T::IndSphereL2) = " ‖⋅‖₂ = r "
fun_par(T::IndSphereL2) = " r = $(round(T.r,3)) "

type NormLinf{T <: AffineOperator} <: NonSmoothTerm
	A::T
	lambda::Real
end

fun_name(T::NormLinf) = " λ ‖⋅‖∞ "
fun_par(T::NormLinf) = " λ = $(round(T.lambda,3)) "

type IndBallLinf{T <: AffineOperator} <: NonSmoothTerm
	A::T
	r::Real
end

fun_name(T::IndBallLinf) = " ‖⋅‖∞ ≤ r "
fun_par(T::IndBallLinf) = " r = $(round(T.r,3)) "

norm(x::OptVar, args...) = norm(eye(x), args...)

function norm{T <: AffineOperator}(A::T, p::Int64=2)
	if p == 0
		return NormL0(A, 1)
	elseif p == 1
		return NormL1(A, 1)
	elseif p == 2
		return NormL2(A, 1)
	else
		error("function not implemented")
	end
end

function norm{T <: AffineOperator}(A::T, p::Float64)
	if p == Inf
		return NormLinf(A, 1)
	else
		error("function not implemented")
	end
end

get_prox(T::NormL0) = ProximalOperators.NormL0(T.lambda)
get_prox(T::NormL1) = ProximalOperators.NormL1(T.lambda)
get_prox(T::NormL2) = ProximalOperators.NormL2(T.lambda)
get_prox(T::NormLinf) = ProximalOperators.NormLinf(T.lambda)

get_prox(T::IndBallL0) = ProximalOperators.IndBallL0(T.r)
get_prox(T::IndBallL1) = ProximalOperators.IndBallL1(T.r)
get_prox(T::IndBallL2) = ProximalOperators.IndBallL2(T.r)
get_prox(T::IndBallLinf) = ProximalOperators.IndBallLinf(T.r)

get_prox(T::IndSphereL2) = ProximalOperators.IndSphereL2(T.r)

*(lambda::Real, T::NormL0) = NormL0(T.A,lambda)
*(lambda::Real,  T::NormL1) = NormL1(T.A,lambda)
*(lambda::Real,  T::NormL2) = NormL2(T.A,lambda)
*(lambda::Real,  T::NormLinf) = NormLinf(T.A,lambda)
 
<=(T::NormL0, r::Integer) = IndBallL0(T.A, r) 
<=(T::NormL1, r::Real)    = IndBallL1(T.A, r/T.lambda) 
<=(T::NormL2, r::Real)    = IndBallL2(T.A, r/T.lambda) 
<=(T::NormLinf, r::Real)    = IndBallLinf(T.A, r/T.lambda) 

==(T::NormL2, r::Real)    = IndSphereL2(T.A, r/T.lambda) 
