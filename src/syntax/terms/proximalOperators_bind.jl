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

# Nuclear norm
function norm(ex::AbstractExpression, ::typeof(*))
	return Term(NuclearNorm(), ex)
end

# Mixed Norm
function norm(ex::AbstractExpression, p1::Int, p2::Int, dim::Int = 1 )
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
		return Term(1.0,Conjugate(Postcompose(t.f,t.lambda)),t.A) 
	else
		error("cannot perform convex conjugation")
	end

end

# Moreau Envelope

export smooth

function smooth(t::Term, gamma = 1.0) 
	return Term(1.0,MoreauEnvelope(Postcompose(t.f,t.lambda),gamma),t.A) 
end

# other stuff, to make Term work with iterators
import Base: start, next, done, isempty
start(t::Term) = false
next(t::Term, state) = (t, true)
done(t::Term, state) =  state
isempty(t::Term) =  false
