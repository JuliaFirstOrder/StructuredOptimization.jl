import Base: size, ~
export Variable

abstract AbstractVariable 

immutable Variable{T <: Union{Real,Complex}} <: AbstractVariable
	x::AbstractArray{T}
end

function Variable(args...)
  Variable(zeros(args...))
end

size(x::Variable) = size(x.x)
"""
`~(x::RegLS.Variable)`

return the `Array` inside `Variable`.
"""
~(x::Variable) = x.x
~{T<:AbstractVariable}(x::Vector{T}) = length(x) == 1 ? ~(x[1]) : (~).(x)

deepcopy!(x::Variable,y::AbstractArray) = deepcopy!(x.x,y) 
deepcopy!{T<:AbstractVariable}(x::Vector{T},y::AbstractArray) = deepcopy!(~x,y) 


