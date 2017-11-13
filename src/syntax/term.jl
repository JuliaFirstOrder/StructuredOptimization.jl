immutable Term{T1 <: Real, T2 <: ProximableFunction, T3 <: AbstractExpression}
	lambda::T1
	f::T2
	A::T3
	Term(lambda::T1, f::T2, ex::T3) where {T1,T2,T3} = new{T1,T2,T3}(lambda,f,ex)
end

function Term{T<:ProximableFunction}(f::T, ex::AbstractExpression)
	A = convert(Expression,ex)
	Term(1,f, A)
end

# Properties

variables(t::Term) = variables(t.A)
operator(t::Term) = operator(t.A)
displacement(t::Term) = displacement(t.A)

#importing properties from AbstractOperators
is_f = [:is_linear,
	:is_eye,
	:is_null,
	:is_diagonal,
	:is_AcA_diagonal,
	:is_AAc_diagonal,
	:is_orthogonal,
	:is_invertible,
	:is_full_row_rank,
	:is_full_column_rank]

for f in is_f
	@eval begin
		import AbstractOperators: $f
		$f(t::Term) = $f(operator(t))
		$f{N}(t::NTuple{N,Term}) = all($f.(t))
	end
end

is_smooth(t::Term) = is_smooth(t.f)

is_convex(t::Term)    = is_convex(t.f) && is_linear(t)

is_quadratic(t::Term) = is_quadratic(t.f) && is_linear(t)

is_strongly_convex(t::Term) = is_strongly_convex(t.f) && is_full_column_rank(operator(t.A))


# Operations

# Define sum of terms simply as their vcat

import Base: +

(+)(a::Term,b::Term) = (a,b)
(+){N}(a::NTuple{N,Term},b::Term) = (a...,b)
(+){N}(a::Term,b::NTuple{N,Term}) = (a,b...)
(+){N}(a::NTuple{N,Term},b::Tuple{}) = a
(+){N}(a::Tuple{},b::NTuple{N,Term}) = b
(+){N,M}(a::NTuple{N,Term},b::NTuple{M,Term}) = (a...,b...)


# Define multiplication by constant

import Base: *

function (*){T1<:Real, T, T2, T3}(a::T1, t::Term{T,T2,T3})
	coeff = *(promote(a,t.lambda)...)
	Term(coeff, t.f, t.A)
end

# Constructors

# Norms

import Base: norm

function norm(ex::AbstractExpression, p=2)
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

# Mixed Norm
function norm(ex::AbstractExpression, p1::Int, p2::Int, dim::Int =1 )
	if p1 == 2 && p2 == 1
		f = NormL21(1.0,dim) 
	else
		error("function not implemented")
	end
	return Term(f, ex)
end

# Norm constraints

import Base: <=, ==

(<=)(t::Term{T1,T2,T3}, r::Integer) where {T1,T2 <: NormL0,T3} =
Term(IndBallL0(round(Int,r/t.lambda)), t.A)
(<=)(t::Term{T1,T2,T3}, r::Real) where {T1, T2 <: NormL1, T3} = Term(IndBallL1(r/t.lambda), t.A)
(<=)(t::Term{T1,T2,T3}, r::Real) where {T1, T2 <: NormL2, T3} = Term(IndBallL2(r/t.lambda), t.A)
(<=)(t::Term{T1,T2,T3}, r::Real) where {T1, T4 <: IndBallL1, T2 <: Conjugate{T4}, T3} = Term(IndBallLinf(r/t.lambda), t.A)

(==)(t::Term{T1,T2,T3}, r::Real)  where {T1,T2 <: NormL2,T3} =    Term(IndSphereL2(r/t.lambda), t.A)

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

(<=)(ex::AbstractExpression, ub) = Term(IndBox(-Inf, ub), ex)
(<=)(lb, ex::AbstractExpression) = Term(IndBox(lb, +Inf), ex)
(>=)(ex::AbstractExpression, lb) = Term(IndBox(lb, +Inf), ex)
(>=)(ub, ex::AbstractExpression) = Term(IndBox(-Inf, ub), ex)

function in(ex::AbstractExpression, bnds::AbstractArray)
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
rank(ex::AbstractExpression) = Term(Rank(), ex)

import Base: <=

(<=)(t::Term{T1,T2,T3} where {T1, T2 <: Rank, T3}, r::Int) = Term(IndBallRank(round(Int,r/t.lambda)), t.A)

# HingeLoss

export hingeloss

hingeloss(ex::AbstractExpression, b::Array{R,1}) where {R <: Real} =
Term(HingeLoss(b), ex)

export sqrhingeloss

sqrhingeloss(ex::AbstractExpression, b::Array{R,1}) where {R <: Real} =
Term(SqrHingeLoss(b), ex)

# CrossEntropy

export crossentropy

crossentropy(ex::AbstractExpression, b::Array{R,1}) where {R <: Real} =
Term(CrossEntropy(b), ex)

# LogBarrier

export logbarrier

logbarrier(ex::AbstractExpression, a::R) where{R <: Real} =
Term(LogBarrier(a), ex)

# HuberLoss

export huberloss

huberloss(ex::AbstractExpression, a::R) where {R <: Real} =
Term(HuberLoss(a), ex)

# Convex conjugate

import Base: conj

function conj(t::Term) 
	if typeof(operator(t)) <: Eye 
		return Term(t.lambda,Conjugate(t.f),t.A) 
	else
		error("cannot perform convex conjugation")
	end

end

# other stuff, to make Term work with iterators
import Base: start, next, done, isempty
start(t::Term) = false
next(t::Term, state) = (t, true)
done(t::Term, state) =  state
isempty(t::Term) =  false
