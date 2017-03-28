import Base: *, <=, ==, isempty
export affine, terms, tilt

immutable CostFunction
	x::Array{OptVar,1}
	f::Array{ExtendedRealValuedFunction,1}
	A::Array{AffineOperator,1}
end

function (cf::CostFunction)(x::AbstractArray, no_affine::Bool = false)
	if no_affine 
		out = 0.0
		for i in eachindex(cf.f)
			out += terms(cf)[i](x) 
		end
		return out
	else
		out = 0.0
		for i in eachindex(cf.f)
			out += terms(cf)[i](affine(cf)[i](x)) 
		end
		return out
	end
end

variable(cf::CostFunction) = cf.x
affine(  cf::CostFunction) = cf.A
terms(   cf::CostFunction) = cf.f
operator(cf::CostFunction) = operator.(cf.A)
tilt(    cf::CostFunction) = tilt.(cf.A)

isempty( cf::CostFunction) = length(cf.f) == 0 ? true : false

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
