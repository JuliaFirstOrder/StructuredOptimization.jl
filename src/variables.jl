import Base:
  size

immutable OptVar{T <: Union{Real,Complex}}
	x::AbstractArray{T}
end

function OptVar(args...)
  OptVar(zeros(args...))
end

size(x::OptVar) = size(x.x)


export OptVar
