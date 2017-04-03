import Base: *, <=, ==, isempty
export affine, terms, tilt

immutable CostFunction
	x::Vector{Variable}
	f::Vector{ExtendedRealValuedFunction}
	A::Vector{AffineOperator}
end

CostFunction() = CostFunction(Vector{Variable}(0),
			      Vector{ExtendedRealValuedFunction}(0),
			      Vector{AffineOperator}(0)) 

variable(cf::CostFunction) = cf.x
affine(  cf::CostFunction) = cf.A
terms(   cf::CostFunction) = cf.f
operator(cf::CostFunction) = operator.(cf.A)
tilt(    cf::CostFunction) = tilt.(cf.A)

isempty( cf::CostFunction) = length(cf.f) == 0 ? true : false
isLeastSquares( cf::CostFunction) =  all([typeof(t) <: LinearLeastSquares for t in terms(cf) ])

function +(h::CostFunction, g::CostFunction) 
	x = addVar(h.x,g.x)
	CostFunction(x,[h.f...,g.f...], [h.A...,g.A...])
end

*(lambda::Real, h::CostFunction)   = CostFunction(h.x,(*).(lambda, h.f),h.A)

<=(h::CostFunction, lambda::Real)  = CostFunction(h.x,(<=).(h.f,lambda),h.A)
==(h::CostFunction, lambda::Real,)  = CostFunction(h.x,(==).(h.f,lambda),h.A)

function cost{T}(cf::CostFunction,resx::Vector{T})
	f = 0.0
	for i in eachindex(terms(cf))
		f += terms(cf)[i](resx[i])
	end
	return f
end

#this function must be used only with sorted and expanded affine operators!
function residual!{T}(resx::Vector{T}, cf::CostFunction, x::AbstractArray)
	f = 0.0
	for i in eachindex(terms(cf))
		evaluate!(resx[i],affine(cf)[i],x)
		f += terms(cf)[i](resx[i])
	end
	return f
end

#this function must be used only with sorted and expanded affine operators!
function residual(cf::CostFunction, x::AbstractArray)
	f = 0.0
	resx = Vector(length(terms(cf)))
	for i in eachindex(terms(cf))
		resx[i] = affine(cf)[i](x)
		f += terms(cf)[i](resx[i])
	end
	return [resx...], f
end

#this function must be used only with sorted and expanded affine operators!
function gradient!{T}(gradfi::Vector{T}, cf::CostFunction, resx::Vector{T} )
	f = 0.0
	for i in eachindex(terms(cf))
		fx, = gradient!(gradfi[i], terms(cf)[i], resx[i])
		f += fx
	end
	return f
end

#this function must be used only with sorted and expanded affine operators!
function gradient{T}(cf::CostFunction, resx::Vector{T} )
	gradfi = deepsimilar(resx)
	f = gradient!(gradfi,cf,resx)
	return gradfi, f
end

#this function must be used only with sorted and expanded affine operators!
#internal variables are used as buffer 
function At_mul_B!{T}(grad::AbstractArray, cf::CostFunction, gradfi::Vector{T} )
	A_mul_B!(grad, adjoint(affine(cf)[1]), gradfi[1]  )
	for i = 2:length(terms(cf))
		A_mul_B!(~variable(cf), adjoint(affine(cf)[i]), gradfi[i] )
		grad .+= ~variable(cf)
	end
end

#this function must be used only with sorted and expanded affine operators!
function At_mul_B{T}(cf::CostFunction, gradfi::Vector{T} )
	grad = deepsimilar(~variable(cf))
	At_mul_B!(grad, cf, gradfi)
	return grad
end

#this function must be used only with sorted and expanded affine operators!
function A_mul_B!{T}(resx::Vector{T}, cf::CostFunction, x::AbstractArray)
	for i in eachindex(terms(cf))
		A_mul_B!(resx[i],operator(affine(cf)[i]),x)
	end
end

function addVar{T}(x::Variable{T},y::Vector)
	any(y.==x) ? y : [y...,x] 
end

function addVar(x::Vector,y::Vector)
	z = copy(y)
	for xi in x
		z = addVar(xi,z)
	end
	return z
end

function Base.show(io::IO, cf::CostFunction)
	if isempty(cf)
		println(io, "Empty Cost Function") 
	else
		description = fun_name(cf.f[1],1)
		operator    = "\n A1 = "*fun_name(RegLS.operator(cf.A[1]))
		parameter   = fun_par(cf.f[1],1)
		for i = 2:length(cf.f)
			description = description*"+ "fun_name(cf.f[i],i)
			operator    = operator*",\n A$i = "*fun_name(RegLS.operator(cf.A[i]))
			parameter = parameter*", "fun_par(cf.f[i],i)
		end
		
		println(io, "description : ", description) 
		println(io, "operators   : ", operator   )
		println(io, "parameters  : ", parameter  )
	end
end
