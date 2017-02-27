import Base: norm, *, <= 

type NormL1{T <: AffineOperator}   <: NonSmoothTerm
	A::T
	lambda::Number
end

type IndBallL1{T <: AffineOperator} <: NonSmoothTerm
	A::T
	r::Real
end

type NormL2{T <: AffineOperator} <: NonSmoothTerm
	A::T
	lambda::Number
end

type IndBallL2{T <: AffineOperator} <: NonSmoothTerm
	A::T
	r::Real
end

norm(x::OptVar, args...) = norm(eye(x), args...)

function norm{T <: AffineOperator}(A::T, p::Int64)
	if p == 1
		return NormL1(A, 1)
	elseif p == 2
		return NormL2(A, 1)
	end
end

function get_prox(T::NormL1)
	return ProximalOperators.NormL1(T.lambda)
end

function get_prox(T::IndBallL1)
	return ProximalOperators.IndBallL1(T.r)
end

function get_prox(T::NormL2)
	return ProximalOperators.NormL2(T.lambda)
end

function get_prox(T::IndBallL2)
	return ProximalOperators.IndBallL2(T.r)
end

*(lambda::Number, T::NormL1) = NormL1(T.A,lambda)
*(lambda::Number, T::NormL2) = NormL2(T.A,lambda)
 
<=(T::NormL1, r::Real) = IndBallL1(T.A, r/T.lambda) 
<=(T::NormL2, r::Real) = IndBallL1(T.A, r/T.lambda) 
