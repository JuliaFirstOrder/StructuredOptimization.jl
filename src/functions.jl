abstract OptTerm
abstract SmoothTerm    <: OptTerm
abstract NonSmoothTerm <: OptTerm

include("functions/get_prox.jl")
include("functions/Norm.jl")
include("functions/Box.jl")
include("functions/LeastSquares.jl")

type CostFunction
	Terms::Array{OptTerm, 1}
end

+(h::OptTerm, g::OptTerm) = CostFunction([h,g])

function +(cf::CostFunction, g::OptTerm) 
	push!(cf.Terms,g)
	return cf
end

