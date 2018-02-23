# utils
export variables, operator, displacement
import Base: convert

convert{T,N,A}(::Type{Expression},x::Variable{T,N,A}) =
Expression{1}((x,),Eye(T,size(x)),zero(T))

"""
`variables(E::Expression)`

Returns a tuple containing the `Variable`s of expression `E`.

# Example

```julia
julia> x,y = Variable(2),Variable(2);

julia> ex = x+y;

julia> variables(ex)
(Variable(Float64, (2,)), Variable(Float64, (2,)))

```

"""
variables(A::Expression)    = A.x
variables(x::Variable)    = x

"""
`operator(E::Expression)`

Returns the `AbstractOperator` of expression `E`.

# Example

```julia
julia> x = Variable(3)
Variable(Float64, (3,))

julia> ex = fft(x);

julia> operator(ex)
ℱ  ℝ^3 -> ℂ^3

```
"""
operator(A::Expression)     = A.L
operator(x::Variable)     = Eye(~x)

"""
`displacement(E::Expression)`

Returns the displacement of expression `E`.

# Example

```julia
julia> x = Variable(3)
Variable(Float64, (3,))

julia> ex = fft(x);

julia> operator(ex)
ℱ  ℝ^3 -> ℂ^3

```
"""
displacement(A::Expression) = A.d
displacement(A::Variable{T}) where {T} = zero(T)
