import Base: *

"""
`*(A::AbstractOperator, ex::AbstractExpression)`

Multiply an 'AbstractExpression` by an `AbstractOperator`. 

# Example

```julia
julia> A = FiniteDiff(Float64, (10,))
δx  ℝ^10 -> ℝ^9

julia> x = Variable(10);

julia> ex = A*x;

julia> B = DCT(9)
ℱc  ℝ^9 -> ℝ^9

julia> ex2 = B*ex;

julia> operator(ex2)
ℱc*δx  ℝ^10 -> ℝ^9
```

"""
function (*)(L::AbstractOperator, a::AbstractExpression)
	A = convert(Expression,a)
	if typeof(displacement(A)) <: Number
		d = displacement(A) == 0. ? zero(codomainType(L)) :
		L*(displacement(A)*ones(codomainType(operator(A)),size(operator(A),1)))
	else
		d = L*displacement(A)
	end
	Expression{length(A.x)}(A.x,L*operator(A),d)
end

"""
`*(A::AbstractMatrix, ex::AbstractExpression)`

Multiply an `AbstractExpression` by an `AbstractMatrix`. 

```julia
julia> A = randn(10,5);

julia> x = Variable(5)
Variable(Float64, (5,))

julia> A*x

```

Other types of multiplications are also possible:

* left array multiplication
  ```julia
  julia> X = Variable(10,5)
  Variable(Float64, (10, 5))
  
  julia> X*randn(5,6)
  
  ```

* scalar multiplication:
  ```julia
  julia> π*X
  
  ```

* elementwise multiplication:
  ```julia
  julia> randn(10,5).*X
  
  ```

"""
function (*)(m::T, a::Union{AbstractVector,AbstractMatrix}) where {T<:AbstractExpression}
	M = convert(Expression,m)
	op = LMatrixOp(codomainType(operator(M)),size(operator(M),1),a)
	return op*M
end
#LMatrixOp

function (*)(M::AbstractMatrix, a::T) where {T<:AbstractExpression}
	A = convert(Expression,a)
	op = MatrixOp(codomainType(operator(A)),size(operator(A),1),M)
	return op*A
end
#MatrixOp

function Base.broadcast(::typeof(*), d::AbstractArray, a::T) where {T<:AbstractExpression}
	A = convert(Expression,a)
	op = DiagOp(codomainType(operator(A)),size(operator(A),1),d)
	return op*A
end
#DiagOp

function (*)(coeff::T1, a::T) where {T1<:Number, T<:AbstractExpression}
	A = convert(Expression,a)
	return Expression{length(A.x)}(A.x,coeff*operator(A),coeff*displacement(A))
end
#Scale

"""
`*(A::AbstractExpression, ex::AbstractExpression)`

Multiply an `AbstractExpression` by another `AbstractExpression`. 

# Examples

```julia
julia> W1 = Variable(10,5)
Variable(Float64, (10, 5))

julia> W2 = Variable(5,15)
Variable(Float64, (5, 15))

julia> ex = W1*σ(W2);

julia> operator(ex)
I*σ  ℝ^(10, 5)  ℝ^(5, 15) -> ℝ^(10, 15)

```

"""
function (*)(ex1::AbstractExpression, ex2::AbstractExpression)
	ex1 = convert(Expression,ex1)
	ex2 = convert(Expression,ex2)
	op = NonLinearCompose(operator(ex1),operator(ex2))
	d = (displacement(ex1) != 0 && displacement(ex2) != 0) ? displacement(ex1)*displacement(ex2) : 
	zero(codomainType(op))
	x = (variables(ex1)...,variables(ex2)...) 
	exp3 = Expression{length(x)}(x,op,d)
	if displacement(ex2) != 0.
		exp3 += Expression{length(ex1.x)}(ex1.x,ex1.L,0.)*displacement(ex2)
	end
	if displacement(ex1) != 0.
		exp3 += displacement(ex1)*Expression{length(ex2.x)}(ex2.x,ex2.L,0.)
	end
	return exp3 
end
# NonLinearCompose
