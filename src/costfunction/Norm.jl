import Base: norm, sum 

norm(x::Variable)    = norm(eye(x))
norm(x::Variable, p) = norm(eye(x), p)

function norm(A::AffineOperator, p::Int64=2)
	if p == 0
		return CostFunction(variable(A), [NormL0(1.)], [A])
	elseif p == 1
		return CostFunction(variable(A), [NormL1(1.)], [A])
	elseif p == 2
		return CostFunction(variable(A), [NormL2(1.)], [A])
	else
		error("function not implemented")
	end
end

function norm(A::AffineOperator, p::Float64)
	if p == Inf
		return CostFunction(variable(A), [NormLinf(1.)], [A])
	else
		error("function not implemented")
	end
end

<=(T::NormL0,   r::Integer) = IndBallL0(r) 
<=(T::NormL1,   r::Real)    = IndBallL1(r/T.lambda) 
<=(T::NormL2,   r::Real)    = IndBallL2(r/T.lambda) 
<=(T::NormLinf, r::Real)  = IndBallLinf(r/T.lambda) 

==(T::NormL2, r::Real)    = IndSphereL2(r/T.lambda) 

sum(T::NormL2,   d::Int)  = NormL21(T.lambda,d) 






