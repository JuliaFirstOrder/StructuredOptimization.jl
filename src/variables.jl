import Base: size
export OptVar, optData

abstract AbstractOptVar 

immutable OptVar{T <: Union{Real,Complex}} <: AbstractOptVar
	x::AbstractArray{T}
end

function OptVar(args...)
  OptVar(zeros(args...))
end

size(x::OptVar) = size(x.x)
optData(x::OptVar) = x.x
optData{T<:AbstractOptVar}(x::Vector{T}) = optData.(x)

deepcopy!(x::OptVar,y::AbstractArray) = deepcopy!(x.x,y) 
deepcopy!{T<:AbstractOptVar}(x::Vector{T},y::AbstractArray) = deepcopy!(optData(x),y) 


