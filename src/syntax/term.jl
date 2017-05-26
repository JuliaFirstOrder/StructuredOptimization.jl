immutable Term{T <: ProximableFunction}
	f::T
	A::AffineExpression
end

Term{T <: ProximableFunction}(f::T, x::Variable) = Term(f, LinearExpression(Eye(eltype(x), size(x)), x))
Term{T <: ProximableFunction}(f::T, L::LinearExpression) = Term{T}(f, AffineExpression([L], zeros(codomainType(L.L), size(L.L, 1))))

# Properties

variables(t::Term) = variables(t.A)
# is_smooth(t::Term) = is_smooth(t.f)
#
# is_smooth(terms::Vararg{Term}) = all(is_smooth.(terms))
#
# is_proximable(t::Term) = length(t.A.Ls) == 1 && is_gram_diagonal(t.A.Ls[1])
#
# is_convex(t::Term) = is_convex(t.f)
#
# is_strongly_convex(t::Term) = is_strongly_convex(t.f) && length(t.A) == 1 && is_full_column_rank(t.A[1])

# Define sum of terms simply as their vcat

import Base: +

(+)(args::Vararg{Term}) = vcat(args...)

# Define multiplication by constant

import Base: *

(*)(a, t::Term) = Term(Postcompose(t.f, a), t.A)

# Norms

import Base: norm

function norm(ex::AbstractAffineExpression, p=2)
	if p == 0
		f = NormL0()
	elseif p == 1
		f = NormL1()
	elseif p == 2
		f = NormL2()
	elseif p == Inf
		f = NormLinf()
	else
		error("function not implemented")
	end
	return Term(f, ex)
end

# Norm constraints

import Base: <=

(<=)(t::Term{T} where T <: NormL0, r::Integer) = Term(IndBallL0(r), t.A)
(<=)(t::Term{T} where T <: NormL1, r::Real) = Term(IndBallL1(r), t.A)
(<=)(t::Term{T} where T <: NormL2, r::Real) = Term(IndBallL2(r), t.A)

# Least square terms

export ls

ls(ex) = Term(SqrNormL2(), ex)

import Base: ^

function (^)(t::Term{T} where T <: NormL2, exp::Integer)
	if exp == 2
		# The coefficient 2.0 is due to the fact that SqrNormL2 divides by 2.0
		return Term(SqrNormL2(2.0), t.A)
	else
		error("function not implemented")
	end
end

# Box constraints

import Base: <=, >=, in

(<=)(ex::AbstractAffineExpression, ub) = Term(IndBox(-Inf, ub), ex)
(<=)(lb, ex::AbstractAffineExpression) = Term(IndBox(lb, +Inf), ex)
(>=)(ex::AbstractAffineExpression, lb) = Term(IndBox(lb, +Inf), ex)
(>=)(ub, ex::AbstractAffineExpression) = Term(IndBox(-Inf, ub), ex)

function in(ex::AbstractAffineExpression, bnds::AbstractArray)
	if length(bnds) != 2
		error("should provide 2 bounds!")
	end
	return Term(IndBox(bnds[1], bnds[2]), ex)
end

# Rank constraints

import Base: rank

# Dirty trick: the "rank" function only makes sense in constraints such as
#   rank(X) <= r,
# therefore here the parameter (1) doesn't really have a role.
# We should probably fix this: it allows weird things in expressing problems.
# Maybe we should have Rank <: ProximableFunction (with no prox! nor gradient!
# defined), that gives IndBallRank when combined with <=.
rank(ex::AbstractAffineExpression) = Term(IndBallRank(1), ex)

import Base: <=

(<=)(t::Term{T} where T <: IndBallRank, r::Integer) = Term(IndBallRank(r), t.A)

# Hinge loss

export hingeloss

hingeloss(x::Variable, args...) = hingeloss(eye(x), args...)

hingeloss{R <: Real}(A::AbstractAffineExpression, b::Array{R,1}) =
Term(variable(A), [HingeLoss(b, 1.0)], [A])
