
# Generalized length, dot product, norm, similar and deepcopy for nested Array objects

deepsimilar{N}(x::NTuple{N,Any}) = similar.(x)

deepsimilar{T <: AbstractArray}(x::T) = similar(x)

deepcopy!{N}(y::NTuple{N,Any}, x::NTuple{N,Any}) = copy.(y,x)

deepcopy!{T <: AbstractArray}(y::T, x::T) = copy!(y,x)

deeplength(x::NTuple) = sum(deeplength.(x))

deeplength{T <: AbstractArray}(x::T) = length(x)

deepvecdot{N}(x::NTuple{N,Any}, y::NTuple{N,Any}) = sum(deepvecdot.(x,y))

deepvecdot{T <: AbstractArray}(x::T, y::T) = vecdot(x, y)

deepvecnorm(x::NTuple) = sqrt(deepvecdot(x, x))

deepvecnorm{T <: AbstractArray}(x::T) = vecnorm(x)

deepmaxabs(x::NTuple) = maximum(abs, deepmaxabs.(x))

deepmaxabs{T <: AbstractArray}(x::T) = maximum(abs,x)

deepmaxabs{T <: Number}(x::T) = abs(x)

deepzeros{N}(t::NTuple{N, DataType}, s::NTuple{N, NTuple}) = zeros.(t, s)

deepzeros{N}(t::Type, s::NTuple{N, Int}) = zeros(t, s)
