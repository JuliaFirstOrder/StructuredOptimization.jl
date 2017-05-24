# norm(x::Variable) = norm(eye(x))
# norm(x::Variable, p) = norm(eye(x), p)
#
# function norm(A::AbstractAffineExpression, p::Int64=2)
# 	if p == 0
# 		return Term(variable(A), [NormL0(1.)], [A])
# 	elseif p == 1
# 		return Term(variable(A), [NormL1(1.)], [A])
# 	elseif p == 2
# 		return Term(variable(A), [NormL2(1.)], [A])
# 	else
# 		error("function not implemented")
# 	end
# end
#
# function norm(A::AbstractAffineExpression, p::Float64)
# 	if p == Inf
# 		return Term(variable(A), [NormLinf(1.)], [A])
# 	else
# 		error("function not implemented")
# 	end
# end

norm(ex::Union{Variable, AbstractAffineExpression}) = norm(ex, 2)

function norm(ex::Union{Variable, AbstractAffineExpression}, p::Integer)
	if p == 0
		f = NormL0()
	elseif p == 1
		f = NormL1()
	elseif p == 2
		f = NormL2()
	else
		error("function not implemented")
	end
	Term(f, ex)
end
