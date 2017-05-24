is_smooth(f::ProximableFunction) = false
is_quadratic(f::ProximableFunction) = false
is_generalized_quadratic(f::ProximableFunction) = false
is_strongly_convex(f::ProximableFunction) = false

include("functions/conjugate.jl")
include("functions/moreauEnvelope.jl")
include("functions/separableSum.jl")
include("functions/sqrNormL2.jl")
include("functions/postcompose.jl")

NullFunction() = IndFree()
is_null(f::ProximableFunction) = false
is_null(f::IndFree) = true

function gradient!(f::ProximableFunction, args...)
	error("gradient not implemented for $f")
end

function gradient(f::ProximableFunction, x::AbstractArray)
	y = similar(x)
	fx = gradient!(y,f,x)
	return y, fx
end

function gradstep!(y::AbstractArray, f::ProximableFunction, x::AbstractArray, gamma)
	fx = gradient!(y, f, x)
	y .*= -gamma
	y .+= x
	return fx
end

function gradstep(f::ProximableFunction, x::AbstractArray, gamma)
	y = similar(x)
	fx = gradstep!(y, f, x, gamma)
	return y, fx
end
