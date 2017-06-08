export gradient!, gradient

is_smooth(f::ProximableFunction) = false
is_quadratic(f::ProximableFunction) = false
is_generalized_quadratic(f::ProximableFunction) = false
is_strongly_convex(f::ProximableFunction) = false

include("functions/conjugate.jl")
include("functions/moreauEnvelope.jl")
include("functions/separableSum.jl")
include("functions/sqrNormL2.jl")
include("functions/PrecomposeDiagonal.jl")
include("functions/postcompose.jl")

function gradient!(f::ProximableFunction, args...)
	error("gradient not implemented for $f")
end

function gradient{T <: Union{AbstractArray, Tuple}}(f::ProximableFunction, x::T)
	y = deepsimilar(x)
	fx = gradient!(y, f, x)
	return y, fx
end

function gradstep!{T <: Union{AbstractArray, Tuple}}(y::T, f::ProximableFunction, x::T, gamma)
	fx = gradient!(y, f, x)
	y .*= -gamma
	y .+= x
	return fx
end

function gradstep{T <: Union{AbstractArray, Tuple}}(f::ProximableFunction, x::T, gamma)
	y = deepsimilar(x)
	fx = gradstep!(y, f, x, gamma)
	return y, fx
end
