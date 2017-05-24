<=(x::Variable, args...) = <=(eye(x), args...)
<={S1 <: Union{Real, AbstractArray}}(lb::S1, x::Variable) = <=(lb, eye(x))

<={S1 <: Union{Real, AbstractArray}}(A::AbstractAffineExpression,  ub::S1) =
Term(variable(A), [IndBox(-Inf, ub)], [A])
<={S1 <: Union{Real, AbstractArray}}(lb::S1,  A::AbstractAffineExpression) =
Term(variable(A), [IndBox(lb, Inf)], [A])

in(x::Variable, args...) = in(eye(x), args...)

function in{S1 <: Union{Real, AbstractArray}}(A::AbstractAffineExpression, bnd::AbstractArray{S1, 1})
	if length(bnd) != 2 error("should provide 2 bounds!") end
	Term(variable(A), [IndBox(bnd[1],bnd[2])], [A])
end
