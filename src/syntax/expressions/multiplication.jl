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

julia> affine(ex2)
ℱc*δx  ℝ^10 -> ℝ^9
```

"""
function (*)(L::AbstractOperator, a::AbstractExpression)
  A = convert(Expression,a)
  Expression{length(A.x)}(A.x,L*affine(A))
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
  op = LMatrixOp(codomainType(affine(M)),size(affine(M),1),a)
  return op*M
end
#LMatrixOp

function (*)(M::AbstractMatrix, a::T) where {T<:AbstractExpression}
  A = convert(Expression,a)
  op = MatrixOp(codomainType(affine(A)),size(affine(A),1),M)
  return op*A
end
#MatrixOp

function Broadcast.broadcasted(::typeof(*), d::D, a::T) where {D <: Union{Number,AbstractArray}, T<:AbstractExpression}
  A = convert(Expression,a)
  op = DiagOp(codomainType(affine(A)),size(affine(A),1),d)
  return op*A
end
Broadcast.broadcasted(::typeof(*), a::T, d::D) where {D <: Union{Number,AbstractArray}, T<:AbstractExpression} = 
d.*a
#DiagOp

function (*)(coeff::T1, a::T) where {T1<:Number, T<:AbstractExpression}
  A = convert(Expression,a)
  return Expression{length(A.x)}(A.x,coeff*affine(A))
end
(*)(a::T, coeff::T1) where {T1<:Number, T<:AbstractExpression} = coeff*a
##Scale

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

julia> affine(ex)
I*σ  ℝ^(10, 5)  ℝ^(5, 15) -> ℝ^(10, 15)

```

`.*(A::AbstractExpression, ex::AbstractExpression)`

Elementwise multiplication between `AbstractExpression` (i.e. Hadamard product). 

"""
function (*)(ex1::AbstractExpression, ex2::AbstractExpression)
  ex1 = convert(Expression,ex1)
  ex2 = convert(Expression,ex2)
  x = extract_variables((ex1,ex2))
  A = extract_affines(x, ex1)
  B = extract_affines(x, ex2)
  op = Ax_mul_Bx(A,B)
  exp3 = Expression{length(x)}(x,op)
  return exp3 
end
# Ax_mul_Bx

function (*)(ex1::AdjointExpression, ex2::AbstractExpression)
  ex1 = ex1.ex
  ex2 = convert(Expression,ex2)
  x = extract_variables((ex1,ex2))
  A = extract_affines(x, ex1)
  B = extract_affines(x, ex2)
  op = Axt_mul_Bx(A,B)
  exp3 = Expression{length(x)}(x,op)
  return exp3 
end
# Axt_mul_Bx

function (*)(ex1::AbstractExpression, ex2::AdjointExpression)
  ex1 = convert(Expression,ex1)
  ex2 = ex2.ex
  x = extract_variables((ex1,ex2))
  A = extract_affines(x, ex1)
  B = extract_affines(x, ex2)
  op = Ax_mul_Bxt(A,B)
  exp3 = Expression{length(x)}(x,op)
  return exp3 
end
# Ax_mul_Bxt

function Broadcast.broadcasted(::typeof(*), ex1::AbstractExpression, ex2::AbstractExpression) 
  ex1 = convert(Expression,ex1)
  ex2 = convert(Expression,ex2)
  x = extract_variables((ex1,ex2))
  A = extract_affines(x, ex1)
  B = extract_affines(x, ex2)
  op = HadamardProd(A,B)
  exp3 = Expression{length(x)}(x,op)
  return exp3 
end
# Hadamard
