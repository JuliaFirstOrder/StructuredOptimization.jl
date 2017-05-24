import Base: size, eltype, ~
export Variable

# TODO: remove for now?
# abstract type AbstractVariable end

immutable Variable{A <: AbstractArray}
	x::A
end

function Variable{I <: Integer}(T::Type, args::Vararg{I})
  Variable(zeros(T, args...))
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
~(x::Tuple{Variable}) = (~).(x)

# deepcopy!(x::Variable,y::AbstractArray) = deepcopy!(x.x,y)
# deepcopy!{T<:AbstractVariable}(x::Vector{T},y::AbstractArray) = deepcopy!(~x,y)

# domainType(x::Variable) = eltype(x.x)

# fun_type(x::AbstractVariable)   =   domainType(x) <: Complex ? "ℂ^$(size(x))" : "ℝ^$(size(x))"
# function fun_type{V <: AbstractVariable}(x::Vector{V})
# 	str = ""
# 	for i in eachindex(x)
# 		str *= fun_type(x[i])
# 		i != length(x) && (str *= ", ")
# 	end
# 	return str
# end
