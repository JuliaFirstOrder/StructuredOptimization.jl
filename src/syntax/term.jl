
immutable Term{T1<:Real, T2 <: ProximableFunction, T3 <: AbstractAffineExpression}
	lambda::T1
	f::T2
	A::T3
	Term{T1,T2,T3}(lambda::T1, f::T2, ex::T3) where {T1,T2,T3} = new{T1,T2,T3}(lambda,f,ex)
end
	
function Term{T<:ProximableFunction}(f::T, ex::AbstractAffineExpression)
	A = convert(AffineExpression,ex)
	Term{Int,T,typeof(A)}(1,f, A)
end


# Properties

variables(t::Term) = variables(t.A)
variables{N}(t::NTuple{N,Term}) = variables.(t)
operator(t::Term) = operator(t.A)
operator{N}(t::NTuple{N,Term}) = operator.(t)
displacement(t::Term) = displacement(t.A)
displacement{N}(t::NTuple{N,Term}) = displacement.(t)

is_smooth(t::Term) = is_smooth(t.f)
is_smooth{N}(t::NTuple{N,Term}) = is_smooth.(t)

is_convex(t::Term) = is_convex(t.f)
is_convex{N}(t::NTuple{N,Term}) = is_convex.(t)

#is_strongly_convex(t::Term) = is_strongly_convex(t.f) && is_full_column_rank(operator(t.A))
#is_strongly_convex{N}(t::NTuple{N,Term}) = is_strongly_convex.(t)

#importing properties from AbstractOperators
is_f = [:is_eye, 
	:is_null, 
	:is_diagonal,
	:is_gram_diagonal, 
	:is_invertible, 
	:is_full_row_rank, 
	:is_full_column_rank]

for f in is_f
	@eval begin
		import AbstractOperators: $f
		$f(t::Term) = $f(operator(t))
		$f{N}(t::NTuple{N,Term}) = $f.(t)
	end
end

# Operations

# Define sum of terms simply as their vcat

import Base: +

(+)(args::Vararg{Term}) = tuple(args...)

# Define multiplication by constant

import Base: *

function (*){T1<:Real, T, T2, T3}(a::T1, t::Term{T,T2,T3})  
	coeff = *(promote(a,t.lambda)...)
	Term{typeof(coeff),T2,T3}(coeff, t.f, t.A)
end

# Constructors

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

(<=)(t::Term{T1,T2,T3} where {T1,T2 <: NormL0,T3}, r::Integer) = 
Term(IndBallL0(round(Int,r/t.lambda)), t.A)
(<=)(t::Term{T1,T2,T3} where {T1,T2 <: NormL1,T3}, r::Real) =    Term(IndBallL1(r/t.lambda), t.A)
(<=)(t::Term{T1,T2,T3} where {T1,T2 <: NormL2,T3}, r::Real) =    Term(IndBallL2(r/t.lambda), t.A)

# Least square terms

export ls

ls(ex) = Term(SqrNormL2(), ex)

import Base: ^

function (^){T1, T2  <: NormL2, T3}(t::Term{T1,T2,T3}, exp::Integer)
	if exp == 2
		# The coefficient 2.0 is due to the fact that SqrNormL2 divides by 2.0
		return t.lambda^2*Term(SqrNormL2(2.0), t.A)
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
immutable Rank <: ProximableFunction end
rank(ex::AbstractAffineExpression) = Term(Rank(), ex)

import Base: <=

(<=)(t::Term{T1,T2,T3} where {T1, T2 <: Rank, T3}, r::Int) = Term(IndBallRank(round(Int,r/t.lambda)), t.A)

# Hinge loss

export hingeloss

hingeloss{R <: Real}(ex::AbstractAffineExpression, b::Array{R,1}) =
Term(HingeLoss(b), ex)



