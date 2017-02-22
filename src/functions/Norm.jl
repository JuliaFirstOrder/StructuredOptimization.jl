import Base: norm, *

type NormL1{T <: AffineOp} <: NonSmoothTerm
	A::T
	lambda::Number
end

norm(x::OptVar, args...) = norm(eye(x), args...)

function norm{T <: AffineOp}(A::T, p::Int64)
	if p == 1
		return NormL1(A, 1)
	end
end

function get_prox(A::NormL1)
	return ProximalOperators.NormL1(A.lambda)
end

*(lambda::Number, T::NormL1) = NormL1(T.A,lambda)
