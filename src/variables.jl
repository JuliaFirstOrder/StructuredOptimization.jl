import Base:
  size

immutable OptVar{T <: Union{Real,Complex}}
	x::AbstractArray{T}
end

function OptVar(args...)
  OptVar(zeros(args...))
end

size(x::OptVar) = size(x.x)
#returns the number of blocks of variables
blkLength(x::OptVar) = 1
blkLength(x::Array{OptVar}) = length(x)


export OptVar
