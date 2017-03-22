import Base: size

abstract AbstractOptVar 

immutable OptVar{T <: Union{Real,Complex}} <: AbstractOptVar
	x::AbstractArray{T}
end

function OptVar(args...)
  OptVar(zeros(args...))
end

size(x::OptVar) = size(x.x)


export OptVar
