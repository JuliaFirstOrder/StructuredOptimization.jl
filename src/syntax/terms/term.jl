struct Term{T1 <: Real, T2 <: ProximableFunction, T3 <: AbstractExpression}
	lambda::T1
	f::T2
	A::T3
	Term(lambda::T1, f::T2, ex::T3) where {T1,T2,T3} = new{T1,T2,T3}(lambda,f,ex)
end

function Term(f::T, ex::AbstractExpression) where {T<:ProximableFunction}
	A = convert(Expression,ex)
	Term(1,f, A)
end

# Operations

# Define sum of terms simply as their vcat

import Base: +

(+)(a::Term,b::Term) = (a,b)
(+)(a::NTuple{N,Term},b::Term)    where {N} = (a...,b)
(+)(a::Term,b::NTuple{N,Term})    where {N} = (a,b...)
(+)(a::NTuple{N,Term},b::Tuple{}) where {N} = a
(+)(a::Tuple{},b::NTuple{N,Term}) where {N} = b
(+)(a::NTuple{N,Term},b::NTuple{M,Term}) where {N,M} = (a...,b...)

# Define multiplication by constant

import Base: *

function (*)(a::T1, t::Term{T,T2,T3}) where {T1<:Real, T, T2, T3}
	coeff = *(promote(a,t.lambda)...)
	Term(coeff, t.f, t.A)
end

function (*)(a::T1, t::T2) where {T1<:Real, N, T2 <: Tuple{Vararg{<:Term,N}} }
	return a.*t 
end

# Properties

variables(t::Term) = variables(t.A)
operator(t::Term) = operator(t.A)
affine(t::Term) = affine(t.A)
displacement(t::Term) = displacement(t.A)

#importing properties from ProximalOperators
import ProximalOperators:
			  is_affine,
			  is_cone,
			  is_convex,
			  is_generalized_quadratic,
			  is_prox_accurate,
			  is_quadratic,
			  is_separable,
			  is_set,
			  is_singleton,
			  is_smooth,
			  is_strongly_convex

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
        :is_full_column_rank,
        :is_sliced
       ]

for f in is_f
	@eval begin
		import AbstractOperators: $f
		$f(t::Term) = $f(operator(t))
		$f(t::NTuple{N,Term}) where {N} = all($f.(t))
	end
end

is_smooth(t::Term) = is_smooth(t.f)
is_convex(t::Term)    = is_convex(t.f) && is_linear(t)
is_quadratic(t::Term) = is_quadratic(t.f) && is_linear(t)
is_strongly_convex(t::Term) = is_strongly_convex(t.f) && is_full_column_rank(operator(t.A))

include("proximalOperators_bind.jl")

# other stuff, to make Term work with iterators
import Base: iterate, isempty
iterate(t::Term, state = true) = state ? (t, false) : nothing
isempty(t::Term) =  false
