import Base: *, <=, ==, isempty
export affine, terms, tilt

immutable CostFunction
	x::Vector{OptVar}
	f::Vector{ExtendedRealValuedFunction}
	A::Vector{AffineOperator}
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

#this function must be used only with sorted and expanded affine operators!
function cost{T}(cf::CostFunction, resx::Vector{T})
	f = 0.0
	for i in eachindex(terms(cf))
		f += terms(cf)[i](resx[i])
	end
	return f
end

#this function must be used only with sorted and expanded affine operators!
function evaluate!{T}(resx::Vector{T}, cf::CostFunction, x::AbstractArray)
	f = 0.0
	for i in eachindex(terms(cf))
		evaluate!(resx[i],affine(cf)[i],x)
		f += terms(cf)[i](resx[i])
	end
	return f
end


#this function must be used only with sorted and expanded affine operators!
function evaluate(cf::CostFunction, x::AbstractArray)
	f = 0.0
	resx = Vector(length(terms(cf)))
	for i in eachindex(terms(cf))
		resx[i] = affine(cf)[i](x)
		f += terms(cf)[i](resx[i])
	end
	return resx, f
end

#this function must be used only with sorted and expanded affine operators!
#internal variables are used as buffer 
function gradient!{T}(grad::AbstractArray, cf::CostFunction, resx::Vector{T} )
	A_mul_B!(grad, adjoint(affine(cf)[1]), resx[1] )
	gradient!(grad, terms(cf)[1])
	for i = 2:length(terms(cf))
		A_mul_B!(optData(variable(cf)), adjoint(affine(cf)[i]), resx[i] )
		gradient!(optData(variable(cf)), terms(cf)[i])
		grad .+= optData(variable(cf))
	end
end

#this function must be used only with sorted and expanded affine operators!
function gradient{T}(cf::CostFunction, resx::Vector{T} )
	grad = deepsimilar(optData(variable(cf)))
	gradient!(grad,cf,resx)
	return grad
end

#this function must be used only with sorted and expanded affine operators!
function shifted_residual!{T}(resx::Vector{T}, cf::CostFunction, x::AbstractArray)
	f = 0.0
	for i in eachindex(terms(cf))
		A_mul_B!(resx[i],operator(affine(cf)[i]),x)
	end
end

function addVar{T}(x::OptVar{T},y::Vector)
	any(y.==x) ? y : [y...,x] 
end

function addVar(x::Vector,y::Vector)
	z = copy(y)
	for xi in x
		z = addVar(xi,z)
	end
	return z
end
