norm(x::Variable) = norm(eye(x))
norm(x::Variable, p) = norm(eye(x), p)

function norm(A::AbstractAffineTerm, p::Int64=2)
	if p == 0
		return CompositeFunction(variable(A), [NormL0(1.)], [A])
	elseif p == 1
		return CompositeFunction(variable(A), [NormL1(1.)], [A])
	elseif p == 2
		return CompositeFunction(variable(A), [NormL2(1.)], [A])
	else
		error("function not implemented")
	end
end

function norm(A::AbstractAffineTerm, p::Float64)
	if p == Inf
		return CompositeFunction(variable(A), [NormLinf(1.)], [A])
	else
		error("function not implemented")
	end
end
