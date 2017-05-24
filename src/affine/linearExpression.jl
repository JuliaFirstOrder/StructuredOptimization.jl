# import Base: +, -
export LinearExpression # operator, adjoint, variable,

immutable LinearExpression{T <: LinearOperator} <: AbstractAffineExpression
	L::T
	x::Variable
	function LinearExpression{T}(L::T, x::Variable) where {T <: LinearOperator}
		if size(L, 2) != size(x)
			throw(ArgumentError("Size of the operator domain $(size(L, 2)) must match size of the variable $(size(x))"))
		end
		if domainType(L) != eltype(x)
			throw(ArgumentError("Type of the operator domain $(domainType(L)) must match type of the variable $(eltype(x))"))
		end
		new(L, x)
	end
end

LinearExpression{T <: LinearOperator}(L::T, x::Variable) = LinearExpression{T}(L, x)

LinearExpression{T <: LinearOperator, E <: LinearExpression}(L::T, e::E) = LinearExpression(Compose(L, e.L), e.x)

variables(A::LinearExpression) = [A.x]

# operator(A::LinearExpression) = A.L
# adjoint(A::LinearExpression) = A.Lt
# tilt(A::LinearExpression) = 0.
#
# domainType(A::LinearExpression) =   domainType(operator(A))
# codomainType(A::LinearExpression) = codomainType(operator(A))
# size(A::LinearExpression, args...) = size(operator(A), args...)

# Constructors

# LinearExpression(x::Variable,L::LinearOperator) = LinearExpression([x],L)

# +(A::LinearExpression) = A
# -(A::LinearExpression) = LinearExpression(variable(A),-operator(A))
#
# +(A::LinearExpression, B::LinearExpression) =  all(variable(A) == variable(B)) ? LinearExpression(variable(A), A.L+B.L) :
# LinearExpression(unsigned_sum(variable(A),operator(A),variable(B),operator(B), true)...)
#
# #see utils down here for unsigned_sum
# # other constructors in AffineConstructors and ComposeAffine
#
# -(A::LinearExpression, B::LinearExpression) =  all(variable(A) == variable(B)) ? LinearExpression(variable(A), A.L-B.L) :
# LinearExpression(unsigned_sum(variable(A),operator(A),variable(B),operator(B), false)...)
#
# -(x::Variable) = -eye(x)
# +(x::Variable,A::LinearExpression) = eye(x)+A
# -(x::Variable,A::LinearExpression) = eye(x)-A
# +(A::LinearExpression,x::Variable) = A+eye(x)
# -(A::LinearExpression,x::Variable) = A-eye(x)
#
# +(x::Variable,y::Variable) = eye(x)+eye(y)
# -(x::Variable,y::Variable) = eye(x)-eye(y)

# Operators

# function evaluate!(y::AbstractArray, A::LinearExpression, x::AbstractArray)
# 	A_mul_B!(y,operator(A),x)
# end
#
# function (A::LinearExpression)(x::AbstractArray)
# 	y = operator(A)*x
# 	return y
# end

#sorting operators

# function sort_and_expand{T<:AbstractVariable}(x::Vector{T}, A::LinearExpression)
# 	if all(x == A.x)
# 		return A
# 	else
# 		dim = size(operator(A),1)
#
# 		H = Vector{LinearOperator}(length(x))
# 		[H[i] = Zeros(dim,size(x[i]))  for i in eachindex(H) ]
# 		for i in eachindex(x)
# 			if any(A.x .== x[i])
# 				idx = find(A.x .== x[i])[1]
# 				H[i] = extract_operator(operator(A),idx)
# 			end
# 		end
# 		H = hcat(H...)
# 		return LinearExpression(x,H)
# 	end
# end

# extract_operator(L::LinearOperator, idx::Int64) = L
# extract_operator(L::HCAT, idx::Int64) = L.A[idx]

# Printing stuff

# fun_name(A::LinearExpression) = fun_name(operator(A))
# fun_type(A::LinearExpression) = fun_type(operator(A))

# utils

# function unsigned_sum(xa::Vector{AbstractVariable}, A::LinearOperator,
# 		      xb::Vector{AbstractVariable}, B::LinearOperator, sign::Bool)
# 	sign ? ([xa[1],xb[1]], hcat(A,B)) : ([xa[1],xb[1]], hcat(A,-B))
# end
#
# function unsigned_sum(xa::Vector{AbstractVariable}, A::HCAT,
# 		      xb::Vector{AbstractVariable}, B::LinearOperator, sign::Bool)
# 	H = copy(A.A)
# 	x = copy(xa)
#
# 	if any(x .== xb[1])
# 		idx = find(x .== xb[1])[1]
# 		sign ? H[idx] = H[idx] + B : H[idx] =  H[idx] - B
# 	else
# 		push!(x,xb[1])
# 		sign ? push!(H,B) : push!(H,-B)
# 	end
# 	return x, HCAT(H,A.mid)
# end
#
# unsigned_sum(xa::Vector{AbstractVariable}, A::LinearOperator,
# 	     xb::Vector{AbstractVariable}, B::HCAT, sign::Bool) = unsigned_sum(xb,B,xa,A,sign)
#
# function unsigned_sum(xa::Vector{AbstractVariable}, A::HCAT,
# 		      xb::Vector{AbstractVariable}, B::HCAT, sign::Bool)
#
# 	x, H = unsigned_sum(xa,A,[xb[1]],B.A[1],sign)
# 	for i = 2:length(xb)
# 		x, H = unsigned_sum(x,H,[xb[i]],B.A[i],sign)
# 	end
# 	return x,H
# end
