# utils
export variables, operator, affine
import Base: convert
import AbstractOperators: displacement

convert(::Type{Expression},x::Variable{T,N,A}) where {T,N,A} =
Expression{1}((x,),Eye(T,size(x)))

"""
`variables(ex::Expression)`

Returns a tuple containing the `Variable`s of expression `ex`.

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
`operator(ex::Expression)`

Returns the `AbstractOperator` of expression `ex`.

# Example

```julia
julia> x = Variable(3)
Variable(Float64, (3,))

julia> ex = fft(x);

julia> operator(ex)
ℱ  ℝ^3 -> ℂ^3

```
"""
operator(A::Expression) = remove_displacement(A.L)
operator(x::Variable)   = Eye(~x)

"""
`affine(ex::Expression)`

Returns the `AbstractOperator` of expression `ex` keeping any affine addition.

"""
affine(A::Expression) = A.L
affine(x::Variable)   = Eye(~x)

"""
`displacement(ex::Expression)`

Returns the displacement of expression `ex`.

# Example

```julia
julia> x = Variable(3)
Variable(Float64, (3,))

julia> ex = fft(x)+[1.+im*2.;0.;3.+im*4];

julia> displacement(ex)
3-element Array{Complex{Float64},1}:
1.0+2.0im
0.0+0.0im
3.0+4.0im

```
"""
displacement(A::Variable{T}) where {T} = zero(T)
displacement(E::Expression) = displacement(affine(E))
