# export affine, terms, tilt, emptycostfun
#
# immutable Term
# 	x::Vector{Variable}
# 	f::Vector{ProximableFunction}
# 	A::Vector{AbstractAffineExpression}
# end
#
# Term() = Term(
# 	Vector{Variable}(0),
#   Vector{ProximableFunction}(0),
#   Vector{AbstractAffineExpression}(0)
# )
#
# variable(cf::Term) = cf.x
# affine(  cf::Term) = cf.A
# terms(   cf::Term) = cf.f
# operator(cf::Term) = operator.(cf.A)
# tilt(    cf::Term) = tilt.(cf.A)
#
# isempty( cf::Term) = length(cf.f) == 0 ? true : false
# # isLeastSquares( cf::Term) =  all([typeof(t) <: LeastSquares for t in terms(cf) ])
#
# function +(h::Term, g::Term)
# 	x = addVar(h.x, g.x)
# 	Term(x, [h.f..., g.f...], [h.A..., g.A...])
# end
#
# *(lambda::Real, h::Term) = Term(h.x, (*).(lambda, h.f), h.A)
#
# function (<=)(h::Term, lambda::Real)
# 	if length(h.f) > 1 error("cannot constrain the sum of functions") end
# 	return Term(h.x, (<=).(h.f, lambda), h.A)
# end
#
# function (==)(h::Term, lambda::Real)
# 	if length(h.f) > 1 error("cannot constrain the sum of functions") end
# 	return Term(h.x, (==).(h.f, lambda), h.A)
# end
#
# sum(h::Term, dim::Int) = Term(h.x, (sum).(h.f, dim), h.A)
#
# function cost{T}(cf::Term,resx::Vector{T})
# 	f = 0.0
# 	for i in eachindex(terms(cf))
# 		f += terms(cf)[i](resx[i])
# 	end
# 	return f
# end
#
# #this function must be used only with sorted and expanded affine operators!
# function residual!{T}(resx::Vector{T}, cf::Term, x::AbstractArray)
# 	# f = 0.0
# 	for i in eachindex(terms(cf))
# 		evaluate!(resx[i],affine(cf)[i],x)
# 		# f += terms(cf)[i](resx[i])
# 	end
# 	# return f
# end
#
# #this function must be used only with sorted and expanded affine operators!
# function residual(cf::Term, x::AbstractArray)
# 	# f = 0.0
# 	resx = Vector(length(terms(cf)))
# 	for i in eachindex(terms(cf))
# 		resx[i] = affine(cf)[i](x)
# 		# f += terms(cf)[i](resx[i])
# 	end
# 	return [resx...] #, f
# end
#
# #this function must be used only with sorted and expanded affine operators!
# function gradient!{T}(gradfi::Vector{T}, cf::Term, resx::Vector{T} )
# 	f = 0.0
# 	for i in eachindex(terms(cf))
# 		fx, = gradient!(gradfi[i], terms(cf)[i], resx[i])
# 		f += fx
# 	end
# 	return f
# end
#
# #this function must be used only with sorted and expanded affine operators!
# function gradient{T}(cf::Term, resx::Vector{T} )
# 	gradfi = deepsimilar(resx)
# 	f = gradient!(gradfi,cf,resx)
# 	return gradfi, f
# end
#
# #this function must be used only with sorted and expanded affine operators!
# #internal variables are used as buffer
# function At_mul_B!{T}(grad::AbstractArray, cf::Term, gradfi::Vector{T} )
# 	A_mul_B!(grad, adjoint(affine(cf)[1]), gradfi[1]  )
# 	for i = 2:length(terms(cf))
# 		A_mul_B!(~variable(cf), adjoint(affine(cf)[i]), gradfi[i] )
# 		grad .+= ~variable(cf)
# 	end
# end
#
# #this function must be used only with sorted and expanded affine operators!
# function At_mul_B{T}(cf::Term, gradfi::Vector{T} )
# 	grad = deepsimilar(~variable(cf))
# 	At_mul_B!(grad, cf, gradfi)
# 	return grad
# end
#
# #this function must be used only with sorted and expanded affine operators!
# function A_mul_B!{T}(resx::Vector{T}, cf::Term, x::AbstractArray)
# 	for i in eachindex(terms(cf))
# 		A_mul_B!(resx[i],operator(affine(cf)[i]),x)
# 	end
# end
#
# function addVar{T}(x::Variable{T},y::Vector)
# 	any(y.==x) ? y : [y...,x]
# end
#
# function addVar(x::Vector,y::Vector)
# 	z = copy(y)
# 	for xi in x
# 		z = addVar(xi,z)
# 	end
# 	return z
# end
#
# function sort_and_expand(x_sorted::Array{Variable,1}, cf::Term)
# 	sA = Vector{AbstractAffineExpression}(length(affine(cf)))
# 	for i in eachindex(affine(cf))
# 		sA[i] = sort_and_expand(x_sorted,affine(cf)[i])
# 	end
# 	return Term(x_sorted,cf.f,sA)
# end
#
# emptycostfun() = Term()
#
# function Base.show(io::IO, cf::Term)
# 	print(io, "Term")
# 	# if isempty(cf)
# 	# 	print(io, "Empty Cost Function")
# 	# else
# 	# 	description = fun_name(cf.f[1],1)
# 	# 	operator    =
# 	# 	"\n A1 = "*fun_name(RegLS.operator(cf.A[1]))*" : "*fun_type(RegLS.operator(cf.A[1]))
# 	# 	parameter   = fun_par(cf.f[1],1)
# 	# 	for i = 2:length(cf.f)
# 	# 		description = description*"+ "fun_name(cf.f[i],i)
# 	# 		operator    *=
# 	# 	"\n A$i = "*fun_name(RegLS.operator(cf.A[i]))*" : "*fun_type(RegLS.operator(cf.A[i]))
# 	# 		parameter = parameter*", "fun_par(cf.f[i],i)
# 	# 	end
# 	#
# 	# 	println(io, "description : ", description)
# 	# 	println(io, "operators   : ", operator   )
# 	# 	print(  io, "parameters  : ", parameter  )
# 	# end
# end
