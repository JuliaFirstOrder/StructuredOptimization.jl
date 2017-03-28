import Base: *, <=, ==

immutable CostFunction
	x::Array{OptVar,1}
	f::Array{ExtendedRealValuedFunction,1}
	A::Array{AffineOperator,1}
end

function +(h::CostFunction, g::CostFunction) 
	x = addVar(h.x,g.x)
	CostFunction(x,[h.f...,g.f...], [h.A...,g.A...])
end

*(lambda::Real, h::CostFunction)   = CostFunction(h.x,(*).(lambda, h.f),h.A)

<=(h::CostFunction, lambda::Real)  = CostFunction(h.x,(<=).(h.f,lambda),h.A)
==(h::CostFunction, lambda::Real,)  = CostFunction(h.x,(==).(h.f,lambda),h.A)


function addVar{T}(x::Vector,y::OptVar{T})
	any(x.==y) ? x : [x...,y] 
end
addVar{T}(x::OptVar{T},y::Vector) = addVar(y,x) 

function addVar(x::Vector,y::Vector)
	z = copy(y)
	for xi in x
		z = addVar(xi,z)
	end
	return z
end
