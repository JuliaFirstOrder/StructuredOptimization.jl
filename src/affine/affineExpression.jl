export AffineExpression

immutable AffineExpression <: AbstractAffineExpression
	Ls::AbstractVector # This will contain LinearExpression objects
	b::AbstractArray
end

# Constructors

AffineExpression(L::LinearExpression) = AffineExpression([L], zeros(codomainType(L.L), size(L.L, 1)))
AffineExpression(Ls::AbstractVector) = AffineExpression(Ls, zeros(codomainType(Ls[1].L), size(Ls[1].L, 1)))
AffineExpression(L::LinearExpression, b::AbstractArray) = AffineExpression([L], b)

variables(A::AffineExpression) = [L.x for L in A.Ls]

# import Base: +, -
#
# immutable AffineExpression <: AbstractAffineExpression
# 	A::LinearExpression
# 	b::AbstractArray
# 	function AffineExpression(A,b)
# 		size(b)!= size(A,1) && throw(DimensionMismatch())
# 		codomainType(A) != eltype(b) && throw(DomainError())
# 		new(A,b)
# 	end
# end
#
# variable(A::AffineExpression) = variable(A.A)
# operator(A::AffineExpression) = operator(A.A)
# adjoint(A::AffineExpression) = adjoint(A.A)
# tilt(A::AffineExpression) = A.b
#
# domainType(A::AffineExpression) = domainType(A::LinearExpression)
# codomainType(A::AffineExpression) = codomainType(A::LinearExpression)
# size(A::AffineExpression, args...) = size(A.A, args...)
#
# # Constructors
#
# (+)(A::LinearExpression,b::AbstractArray) = AffineExpression( A, b)
# (-)(A::LinearExpression,b::AbstractArray) = AffineExpression( A,-b)
# (+)(b::AbstractArray,A::LinearExpression) = AffineExpression( A, b)
# (-)(b::AbstractArray,A::LinearExpression) = AffineExpression(-A, b)
# (+)(A::AffineExpression,b::AbstractArray) = AffineExpression( A.A, A.b+b)
# (-)(A::AffineExpression,b::AbstractArray) = AffineExpression( A.A, A.b-b)
# (+)(b::AbstractArray,A::AffineExpression) = AffineExpression( A.A, A.b+b)
# (-)(b::AbstractArray,A::AffineExpression) = AffineExpression(-A.A, b-A.b)
#
# (+)(A::AffineExpression,B::LinearExpression) = AffineExpression(A.A+B, A.b )
# (-)(A::AffineExpression,B::LinearExpression) = AffineExpression(A.A-B, A.b )
# (+)(A::LinearExpression,B::AffineExpression) = AffineExpression(A+B.A, B.b )
# (-)(A::LinearExpression,B::AffineExpression) = AffineExpression(A-B.A,-B.b )
#
# (+)(A::AffineExpression,B::AffineExpression) = AffineExpression(A.A+B.A, A.b+B.b )
# (-)(A::AffineExpression,B::AffineExpression) = AffineExpression(A.A-B.A, A.b-B.b )
#
# (+)(x::Variable,b::AbstractArray) = eye(x)+b
# (-)(x::Variable,b::AbstractArray) = eye(x)-b
# (+)(b::AbstractArray,x::Variable) = b+eye(x)
# (-)(b::AbstractArray,x::Variable) = b-eye(x)
#
# #special cases
#
# (+)(x::Variable,b::Number) = b == 0 ? eye(x) : error("cannot sum $(typeof(x)) with $(typeof(b))")
# (-)(x::Variable,b::Number) = b == 0 ? eye(x) : error("cannot sum $(typeof(x)) with $(typeof(b))")
# (+)(A::AbstractAffineExpression,b::Number) = b == 0 ? A : error("cannot sum $(typeof(A)) with $(typeof(b))")
# (-)(A::AbstractAffineExpression,b::Number) = b == 0 ? A : error("cannot sum $(typeof(A)) with $(typeof(b))")
#
# # Operators
#
# function evaluate!(y::AbstractArray, A::AffineExpression, x::AbstractArray)
# 	evaluate!(y, A.A, x)
# 	y .+= A.b
# end
#
# function (A::AffineExpression)(x::AbstractArray)
# 	y = A.A(x)
# 	y .+= A.b
# 	return y
# end
#
# #sorting operators
#
# sort_and_expand(x, A::AffineExpression) = AffineExpression(sort_and_expand(x, A.A),A.b)
