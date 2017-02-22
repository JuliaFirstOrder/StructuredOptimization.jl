abstract OptTerm
abstract SmoothTerm    <: OptTerm
abstract NonSmoothTerm <: OptTerm

import Base: norm, *

type NormL1{T<:AffineOp} <: NonSmoothTerm
	A::T
	lambda::Number
end

function norm(x::OptVar, p::Int64)
	if p == 1
		return NormL1(eye(x),1) 
	end
end

function norm{T<:AffineOp}(A::T, p::Int64)
	if p == 1
		return NormL1(A,1) 
	end
end

function get_prox(A::NormL1)
	return ProximalOperators.NormL1(A.lambda)
end

type LeastSquares{T<:AffineOp} <: SmoothTerm
	A::T
end

ls{T<:AffineOp}(A::T) = LeastSquares(A)

export ls


*(lambda::Number, T::NormL1) = NormL1(T.A,lambda)  

type CostFunction
	Terms::Array{OptTerm,1}
end

+(h::OptTerm,g::OptTerm) = CostFunction([h,g])
+(cf::CostFunction, g::OptTerm) = push!(cf.Terms,g)




