import Base: convert, size, eltype, ~
export Variable

immutable Variable{T, N, A <: AbstractArray{T,N}} <: AbstractExpression
	x::A
end

# constructors
"""
`Variable([T::Type,] dims...)`

Returns a `Variable` of dimension `dims` initialized with an array of all zeros.

`Variable(x::AbstractArray)`

Returns a `Variable` of dimension `size(x)` initialized with `x`

"""
Variable{T,N,A<:AbstractArray{T,N}}(x::A) = Variable{T,N,A}(x)

function Variable{I <: Integer,N}(T::Type, args::Vararg{I,N})
	Variable{T,N,Array{T,N}}(zeros(T, args...))
end

function Variable{I <: Integer}(args::Vararg{I})
  Variable(zeros(args...))
end

# Utils

function Base.show(io::IO, x::Variable)
  print(io, "Variable($(eltype(x.x)), $(size(x.x)))")
end


"""
`~(x::Variable)`

Returns the `Array` of the variable `x`
"""
~(x::Variable) = x.x
~(x::Tuple{Variable}) = (~)(x[1])
~{N}(x::NTuple{N,Variable}) = (~).(x)

"""
size(x::Variable, [dim...])

Like `size(A::AbstractArray, [dims...])` returns the tuple containing the dimensions of the variable `x`.
"""
size(x::Variable) = size(x.x)
size(x::Variable, dim::I) where { I <: Integer} = size(x.x, dim)

"""
eltype(x::Variable)

Like `eltype(x::AbstractArray)` returns the type of the elements of `x`.
"""
eltype(x::Variable) = blockeltype(x.x)
