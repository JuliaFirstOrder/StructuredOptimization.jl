export AffineExpression

immutable AffineExpression <: AbstractAffineExpression
	Ls::Vector{LinearExpression} # This will contain LinearExpression objects
	b::AbstractArray
	function AffineExpression(myLs, myb)
		if size(myLs[1].L, 1) != size(myb)
			error("size of b must match that of linear terms")
		end
		if codomainType(myLs[1].L) != eltype(myb)
			error("type of b must match that of linear terms")
		end
		if length(myLs) > 1 && any([size(myLs[1].L, 1) != size(myLs[k].L, 1) for k in 2:length(myLs)])
			error("size of all linear terms must be the same")
		end
		if length(myLs) > 1 && any([codomainType(myLs[1].L) != codomainType(myLs[k].L) for k in 2:length(myLs)])
			error("type of all linear terms must be the same")
		end
		new(myLs, myb)
	end
end

# Constructors

AffineExpression(L::LinearExpression) = AffineExpression([L], zeros(codomainType(L.L), size(L.L, 1)))
AffineExpression(Ls::AbstractVector) = AffineExpression(Ls, zeros(codomainType(Ls[1].L), size(Ls[1].L, 1)))
AffineExpression(L::LinearExpression, b::AbstractArray) = AffineExpression([L], b)

# Properties

variables(A::AffineExpression) = Set([L.x for L in A.Ls])

# Variable unary operations

import Base: +, -

(+)(x::Variable) = LinearExpression(Eye(eltype(x), size(x)), x)
(-)(x::Variable) = LinearExpression(-Eye(eltype(x), size(x)), x)

# Variable +/- AbstractArray

(+)(x::Variable, b::AbstractArray) = AffineExpression([LinearExpression(Eye(eltype(x), size(x)), x)], b)
(-)(x::Variable, b::AbstractArray) = AffineExpression([LinearExpression(Eye(eltype(x), size(x)), x)], -b)
(+)(b::AbstractArray, x::Variable) = (+)(x, b)
(-)(b::AbstractArray, x::Variable) = AffineExpression([LinearExpression(-Eye(eltype(x), size(x)), x)], -b)

# Variable +/- Variable

(+)(x::Variable, y::Variable) = AffineExpression([LinearExpression(Eye(eltype(x), size(x)), x), LinearExpression(Eye(eltype(y), size(y)), y)])
(-)(x::Variable, y::Variable) = AffineExpression([LinearExpression(Eye(eltype(x), size(x)), x), LinearExpression(-Eye(eltype(y), size(y)), y)])

# LinearExpression unary operations

(+)(ex::LinearExpression) = ex1
(-)(ex::LinearExpression) = LinearExpression(-ex.L, ex.x)

# LinearExpression +/- AbstractArray

(+)(ex::LinearExpression, b::AbstractArray) = AffineExpression([ex], b)
(-)(ex::LinearExpression, b::AbstractArray) = AffineExpression([ex], -b)
(+)(b::AbstractArray, ex::LinearExpression) = (+)(ex, b)
(-)(b::AbstractArray, ex::LinearExpression) = AffineExpression([-ex], b)

# LinearExpression +/- Variable

(+)(ex::LinearExpression, y::Variable) = AffineExpression([ex, LinearExpression(Eye(eltype(y), size(y)), y)])
(-)(ex::LinearExpression, y::Variable) = AffineExpression([ex, LinearExpression(-Eye(eltype(y), size(y)), y)])
(+)(y::Variable, ex::LinearExpression) = (+)(ex, y)
(-)(y::Variable, ex::LinearExpression) = AffineExpression([LinearExpression(Eye(eltype(y), size(y)), y), -ex])

# LinearExpression +/- LinearExpression

(+)(lex1::LinearExpression, lex2::LinearExpression) = AffineExpression([lex1, lex2])
(-)(lex1::LinearExpression, lex2::LinearExpression) = AffineExpression([lex1, -lex2])

# AffineExpression +/- AbstractArray

(+)(ex::AffineExpression, b::AbstractArray) = AffineExpression(ex.Ls, ex.b+b)
(-)(ex::AffineExpression, b::AbstractArray) = AffineExpression(ex.Ls, ex.b-b)
(+)(b::AbstractArray, ex::AffineExpression) = (+)(ex, b)
(-)(b::AbstractArray, ex::AffineExpression) = AffineExpression(.-ex.Ls, b-ex.b)

# AffineExpression +/- Variable

(+)(ex::AffineExpression, y::Variable) = AffineExpression(vcat(ex.Ls, LinearExpression(Eye(eltype(y), size(y)), y)), ex.b)
(-)(ex::AffineExpression, y::Variable) = AffineExpression(vcat(ex.Ls, LinearExpression(-Eye(eltype(y), size(y)), y)), ex.b)
(+)(y::Variable, ex::AffineExpression) = (+)(ex, y)
(-)(y::Variable, ex::AffineExpression) = AffineExpression(vcat(.-ex.Ls, LinearExpression(Eye(eltype(y), size(y)), y)), b)

# AffineExpression +/- LinearExpression

(+)(ex::AffineExpression, lex::LinearExpression) = AffineExpression(vcat(ex.Ls, lex), ex.b)
(-)(ex::AffineExpression, lex::LinearExpression) = AffineExpression(vcat(ex.Ls, -lex), ex.b)
(+)(lex::LinearExpression, ex::AffineExpression) = (+)(ex, lex)
(-)(lex::LinearExpression, ex::AffineExpression) = AffineExpression(vcat(.-ex.Ls, lex), ex.b)

# AffineExpression +/- AffineExpression

(+)(ex1::AffineExpression, ex2::AffineExpression) = AffineExpression(vcat(ex1.Ls, ex2.Ls), ex1.b+ex2.b)
(-)(ex1::AffineExpression, ex2::AffineExpression) = AffineExpression(vcat(ex1.Ls, .-ex2.Ls), ex1.b-ex2.b)
