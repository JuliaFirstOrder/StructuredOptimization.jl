import Base: convert, size, eltype, ndims, ~
export Variable

immutable Variable{T, N, A <: AbstractArray{T,N}} <: AbstractAffineExpression
	x::A
end

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

size(x::Variable) = size(x.x)
size(x::Variable, i) = size(x.x, i)
eltype(x::Variable) = eltype(x.x)
ndims(x::Variable) = ndims(x.x)

"""
  `~(x::RegLS.Variable)`

returns the `Array` object containing the value of `x`.
"""
~(x::Variable) = x.x
~(x::Tuple{Variable}) = (~)(x[1])
~{N}(x::NTuple{N,Variable}) = (~).(x)


# other stuff, to make Variable work with iterators
import Base: start, next, done, isempty
start(t::Variable) = false
next(t::Variable, state) = (t, true)
done(t::Variable, state) =  state
isempty(t::Variable) =  false


