# Norms

import LinearAlgebra: norm
export norm

"""
`norm(x::AbstractExpression, p=2, [q,] [dim=1])`

Returns the norm of `x`. 

Supported norms:

* `p = 0` ``l_0``-pseudo-norm

* `p = 1` ``l_1``-norm

* `p = 2` ``l_2``-norm

* `p = Inf` ``l_{\\infty}``-norm

* `p = *` nuclear norm

* `p = 2`, `q = 1` ``l_{2,1}`` mixed norm (aka Sum-of-``l_2``-norms) 
```math
f(\\mathbf{X}) = \\sum_i \\| \\mathbf{x}_i \\|
```
where ``\\mathbf{x}_i`` is the ``i``-th column if `dim == 1` (or row if  `dim == 2`) of ``\\mathbf{X}``.

"""
function norm(ex::AbstractExpression, p::Real=2)
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

# Least square terms

export ls

"""
`ls(x::AbstractExpression)`

Returns the squared norm (least squares) of `x`:

```math
f (\\mathbf{x}) = \\frac{1}{2} \\| \\mathbf{x} \\|^2
```

(shorthand of `1/2*norm(x)^2`).
"""
ls(ex) = Term(SqrNormL2(), ex)

import Base: ^

function (^)(t::Term{T1,T2,T3}, exp::Integer) where {T1, T2  <: NormL2, T3}
	if exp == 2
		# The coefficient 2.0 is due to the fact that SqrNormL2 divides by 2.0
		return t.lambda^2*Term(SqrNormL2(2.0), t.A)
	else
		error("function not implemented")
	end
end

# HingeLoss

export hingeloss

"""
`hingeloss(x::AbstractExpression, y::Array)`

Applies the Hinge loss function 
```math
f( \\mathbf{x} ) = \\sum_{i} \\max\\{0, 1 - y_i x_i \\},
```
where `y` is an array containing ``y_i``.
"""
hingeloss(ex::AbstractExpression, b::Array{R,1}) where {R <: Real} =
Term(HingeLoss(b), ex)

# HingeLoss

export sqrhingeloss

"""
`sqrhingeloss(x::AbstractExpression, y::Array)`

Applies the squared Hinge loss function 
```math
f( \\mathbf{x} ) = \\sum_{i} \\max\\{0, 1 - y_i x_i \\}^2,
```
where `y` is an array containing ``y_i``.
"""
sqrhingeloss(ex::AbstractExpression, b::Array{R,1}) where {R <: Real} =
Term(SqrHingeLoss(b), ex)

# CrossEntropy

export crossentropy

"""
`crossentropy(x::AbstractExpression, y::Array)`

Applies the cross entropy loss function: 
```math
f(\\mathbf{x}) = -1/N \\sum_{i}^{N} y_i \\log (x_i)+(1-y_i) \\log (1-x_i),
```
where `y` is an array of length ``N`` containing ``y_i`` having ``0 \\leq y_i \\leq 1``.
"""
crossentropy(ex::AbstractExpression, b::Array{R,1}) where {R <: Real} =
Term(CrossEntropy(b), ex)

# LogisticLoss

export logisticloss

"""
`logbarrier(x::AbstractExpression, y::AbstractArray)`

Applies the logistic loss function: 
```math
f(\\mathbf{x}) = \\sum_{i} \\log(1+ \\exp(-y_i x_i)), 
```
where `y` is an array containing ``y_i``.
"""
logisticloss(ex::AbstractExpression, y::AbstractArray) =
Term(LogisticLoss(y, 1.0), ex)

# LogBarrier

export logbarrier

"""
`logbarrier(x::AbstractExpression)`

Applies the log barrier function: 
```math
f(\\mathbf{x}) = -\\sum_i \\log( x_i ).
```
"""
logbarrier(ex::AbstractExpression) =
Term(LogBarrier(1.0), ex)

# HuberLoss

export huberloss

"""
`huberloss(x::AbstractExpression, ρ=1.0)`

Applies the Huber loss function: 
```math
f(\\mathbf{x}) = \\begin{cases}
  \\tfrac{1}{2}\\| \\mathbf{x} \\|^2 & \\text{if} \\ \\| \\mathbf{x} \\| \\leq \\rho \\\\
  \\rho (\\| \\mathbf{x} \\| - \\tfrac{\\rho}{2}) & \\text{otherwise}.
\\end{cases}
```
"""
huberloss(ex::AbstractExpression, rho::R = 1.0) where {R <: Real} =
Term(HuberLoss(rho), ex)

import Base: maximum

"""
`maximum(x::AbstractExpression)`

Applies the function: 
```math
f(\\mathbf{x}) = \\max \\{x_i : i = 1,\\ldots, n \\}.
```
"""
maximum(ex::AbstractExpression) =
Term(Maximum(), ex)

export sumpositive

"""
`sumpositive(x::AbstractExpression, ρ=1.0)`

Applies the function: 
```math
f(\\mathbf{x}) = \\sum_i \\max \\{x_i, 0\\}.
```
"""
sumpositive(ex::AbstractExpression) =
Term(SumPositive(), ex)

import LinearAlgebra: dot
export dot

"""
`dot(c::AbstractVector, x::AbstractExpression)`

Applies the function: 
```math
f(\\mathbf{x}) = \\mathbf{c}^{T}\\mathbf{x}.
```
"""
dot(c::AbstractVector, ex::AbstractExpression) =
Term(Linear(c), ex)


# Inequalities

import Base: <=
"""
Inequalities constrains 

## Norm Inequalities constraints

* `norm(x::AbstractExpression, 0) <= n::Integer` 

  ``\\mathrm{nnz}(\\mathbf{x}) \\leq n``

* `norm(x::AbstractExpression, 1) <= r::Number` 

  ``\\sum_i \\| x_i \\| \\leq r``

* `norm(x::AbstractExpression, 2) <= r::Number` 

  ``\\| \\mathbf{x} \\| \\leq r``

* `norm(x::AbstractExpression, Inf) <= r::Number` 

  `` \\max \\{ x_1, x_2, \\dots \\}  \\leq r``

## Box inequality constraints

* `x::AbstractExpression <= u::Union{AbstractArray, Real}`

  `` x_i \\leq u_i ``

* `x::AbstractExpression >= l::Union{AbstractArray, Real}`

  `` x_i \\geq l_i ``

  Notice that the notation `x in [l,u]` is also possible.

## Rank inequality constraints

* `rank(X::AbstractExpression) <= n::Integer`

  ``\\mathrm{rank}(\\mathbf{X}) \\leq r`` 

  Notice that the expression `X` must have a codomain with dimension equal to 2. 

"""
(<=)(t::Term{T1,T2,T3}, r::Integer) where {T1,T2 <: NormL0,T3} =
Term(IndBallL0(round(Int,r/t.lambda)), t.A)
(<=)(t::Term{T1,T2,T3}, r::Real) where {T1, T2 <: NormL1, T3} = Term(IndBallL1(r/t.lambda), t.A)
(<=)(t::Term{T1,T2,T3}, r::Real) where {T1, T2 <: NormL2, T3} = Term(IndBallL2(r/t.lambda), t.A)
(<=)(t::Term{T1,T2,T3}, r::Real) where {T1, T4 <: IndBallL1, T2 <: Conjugate{T4}, T3} = Term(IndBallLinf(r/t.lambda), t.A)

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

import LinearAlgebra: rank
export rank

# Dirty trick: the "rank" function only makes sense in constraints such as
#   rank(X) <= r,
# therefore here the parameter (1) doesn't really have a role.
# We should probably fix this: it allows weird things in expressing problems.
# Maybe we should have Rank <: ProximableFunction (with no prox! nor gradient!
# defined), that gives IndBallRank when combined with <=.
struct Rank <: ProximableFunction end
rank(ex::AbstractExpression) = Term(Rank(), ex)

import Base: <=

(<=)(t::Term{T1,T2,T3} where {T1, T2 <: Rank, T3}, r::Int) = Term(IndBallRank(round(Int,r/t.lambda)), t.A)

import Base: ==

"""
Equalities constraints

## Affine space constraint

* `ex == b::Union{Real,AbstractArray}` 

  Requires expression to be affine.

  ### Example
  ```julia
  julia> A,b  = randn(10,5), randn(10);

  julia> x = Variable(5);
  
  julia> A*x == b
  ```

## Norm equality constraint

* `norm(x::AbstractExpression) == r::Number` 

  ``\\| \\mathbf{x} \\| = r``

## Binary constraint

* `x::AbstractExpression == (l, u)` 

  ``\\mathbf{x} = \\mathbf{l}`` or ``\\mathbf{x} = \\mathbf{u}``

"""
(==)(t::Term{T1,T2,T3}, r::Real)  where {T1,T2 <: NormL2,T3} = Term(IndSphereL2(r/t.lambda), t.A)
# IndSphereL2

(==)(ex::AbstractExpression, lu::Tuple{Union{Real,AbstractArray},Union{Real,AbstractArray}}) = 
Term(IndBinary(lu...), ex)
# IndBinary

# IndAffine
function (==)(ex::AbstractExpression, b::Union{Real,AbstractArray}) 
    op = operator(ex)
    d  = displacement(ex)
    if typeof(op) <: MatrixOp
        A = op.A
        bb = b.-d
        p = IndAffine(A, bb)
        return Term(p, variables(ex)[1])
    else
       # TODO change this
       error("Currently affine equality supported only with `MatrixOp`")
    end
end

# Transforms
# Convex conjugate
import Base: conj

"""
`conj(t::Term)`

Returns the convex conjugate transform of `t`:
```math
f^*(\\mathbf{x}) = \\sup_{\\mathbf{y}} \\{ \\langle \\mathbf{y}, \\mathbf{x} \\rangle - f(\\mathbf{y}) \\}.
```

# Example
```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> x = Variable(4);

julia> t = conj(norm(x,1))

```

"""
function conj(t::Term) 
	if typeof(operator(t)) <: Eye 
		return Term(1.0,Conjugate(Postcompose(t.f,t.lambda)),t.A) 
	else
		error("cannot perform convex conjugation")
	end

end

# Moreau Envelope
export smooth

"""
`smooth(t::Term, gamma = 1.0)`

Smooths the nonsmooth term `t` using Moreau envelope:

```math
f^{\\gamma}(\\mathbf{x}) = \\min_{\\mathbf{z}} \\left\\{ f(\\mathbf{z}) + \\tfrac{1}{2\\gamma}\\|\\mathbf{z}-\\mathbf{x}\\|^2 \\right\\}.
```

# Example
```julia
julia> x = Variable(4)
Variable(Float64, (4,))

julia> x = Variable(4);

julia> t = smooth(norm(x,1))

```

"""
function smooth(t::Term, gamma = 1.0) 
    if !is_smooth(t)
        return Term(1.0,MoreauEnvelope(Postcompose(t.f,t.lambda),gamma),t.A) 
    else
        return t
    end
end
