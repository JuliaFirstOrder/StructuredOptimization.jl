import Base: *, <=, ==, sum, isempty
export affine, terms, tilt

immutable CompositeFunction
	x::Vector{Variable}
	f::Vector{ExtendedRealValuedFunction}
	A::Vector{AbstractAffineTerm}
end

CompositeFunction() = CompositeFunction(Vector{Variable}(0),
			      Vector{ExtendedRealValuedFunction}(0),
			      Vector{AbstractAffineTerm}(0))

variable(cf::CompositeFunction) = cf.x
affine(  cf::CompositeFunction) = cf.A
terms(   cf::CompositeFunction) = cf.f
operator(cf::CompositeFunction) = operator.(cf.A)
tilt(    cf::CompositeFunction) = tilt.(cf.A)

isempty( cf::CompositeFunction) = length(cf.f) == 0 ? true : false
isLeastSquares( cf::CompositeFunction) =  all([typeof(t) <: LeastSquares for t in terms(cf) ])

function +(h::CompositeFunction, g::CompositeFunction)
	x = addVar(h.x,g.x)
	CompositeFunction(x,[h.f...,g.f...], [h.A...,g.A...])
end

*(lambda::Real, h::CompositeFunction)   = CompositeFunction(h.x,(*).(lambda, h.f),h.A)

<=(h::CompositeFunction, lambda::Real)  = CompositeFunction(h.x,(<=).(h.f,lambda),h.A)
==(h::CompositeFunction, lambda::Real)  = CompositeFunction(h.x,(==).(h.f,lambda),h.A)
sum(h::CompositeFunction,  dim::Int  ) = CompositeFunction(h.x,(sum).(h.f,dim   ),h.A)

function cost{T}(cf::CompositeFunction,resx::Vector{T})
	f = 0.0
	for i in eachindex(terms(cf))
		f += terms(cf)[i](resx[i])
	end
	return f
end

#this function must be used only with sorted and expanded affine operators!
function residual!{T}(resx::Vector{T}, cf::CompositeFunction, x::AbstractArray)
	f = 0.0
	for i in eachindex(terms(cf))
		evaluate!(resx[i],affine(cf)[i],x)
		f += terms(cf)[i](resx[i])
	end
	return f
end

#this function must be used only with sorted and expanded affine operators!
function residual(cf::CompositeFunction, x::AbstractArray)
	f = 0.0
	resx = Vector(length(terms(cf)))
	for i in eachindex(terms(cf))
		resx[i] = affine(cf)[i](x)
		f += terms(cf)[i](resx[i])
	end
	return [resx...], f
end

#this function must be used only with sorted and expanded affine operators!
function gradient!{T}(gradfi::Vector{T}, cf::CompositeFunction, resx::Vector{T} )
	f = 0.0
	for i in eachindex(terms(cf))
		fx, = gradient!(gradfi[i], terms(cf)[i], resx[i])
		f += fx
	end
	return f
end

#this function must be used only with sorted and expanded affine operators!
function gradient{T}(cf::CompositeFunction, resx::Vector{T} )
	gradfi = deepsimilar(resx)
	f = gradient!(gradfi,cf,resx)
	return gradfi, f
end

#this function must be used only with sorted and expanded affine operators!
#internal variables are used as buffer
function At_mul_B!{T}(grad::AbstractArray, cf::CompositeFunction, gradfi::Vector{T} )
	A_mul_B!(grad, adjoint(affine(cf)[1]), gradfi[1]  )
	for i = 2:length(terms(cf))
		A_mul_B!(~variable(cf), adjoint(affine(cf)[i]), gradfi[i] )
		grad .+= ~variable(cf)
	end
end

#this function must be used only with sorted and expanded affine operators!
function At_mul_B{T}(cf::CompositeFunction, gradfi::Vector{T} )
	grad = deepsimilar(~variable(cf))
	At_mul_B!(grad, cf, gradfi)
	return grad
end

#this function must be used only with sorted and expanded affine operators!
function A_mul_B!{T}(resx::Vector{T}, cf::CompositeFunction, x::AbstractArray)
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

function sort_and_expand(x_sorted::Array{Variable,1}, cf::CompositeFunction)
	sA = Vector{AbstractAffineTerm}(length(affine(cf)))
	for i in eachindex(affine(cf))
		sA[i] = sort_and_expand(x_sorted,affine(cf)[i])
	end
	return CompositeFunction(x_sorted,cf.f,sA)
end

function Base.show(io::IO, cf::CompositeFunction)
	if isempty(cf)
		print(io, "Empty Cost Function")
	else
		description = fun_name(cf.f[1],1)
		operator    =
		"\n A1 = "*fun_name(RegLS.operator(cf.A[1]))*" : "*fun_type(RegLS.operator(cf.A[1]))
		parameter   = fun_par(cf.f[1],1)
		for i = 2:length(cf.f)
			description = description*"+ "fun_name(cf.f[i],i)
			operator    *=
		"\n A$i = "*fun_name(RegLS.operator(cf.A[i]))*" : "*fun_type(RegLS.operator(cf.A[i]))
			parameter = parameter*", "fun_par(cf.f[i],i)
		end

		println(io, "description : ", description)
		println(io, "operators   : ", operator   )
		print(  io, "parameters  : ", parameter  )
	end
end
