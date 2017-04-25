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

domainType(x::Variable) = eltype(~x) 


function Base.show{V <: AbstractVariable }(io::IO, f::V)
  println(io, "description : Optization Variable")
  println(io, "type        : ", fun_type(f))
end

fun_type(x::AbstractVariable)   =   domainType(x) <: Complex ? "ℂ^$(size(x))" : "ℝ^$(size(x))"
function fun_type{V <: AbstractVariable}(x::Vector{V})     
	str = ""
	for i in eachindex(x) 
		str *= fun_type(x[i]) 
		i != length(x) && (str *= ", ")
	end
	return str
end

