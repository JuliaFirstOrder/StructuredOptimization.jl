abstract OptTerm
abstract SmoothTerm    <: OptTerm
abstract NonSmoothTerm <: OptTerm

include("functions/Norm.jl")
include("functions/LeastSquares.jl")

type CostFunction
	Terms::Array{OptTerm, 1}
end

+(h::OptTerm, g::OptTerm) = CostFunction([h,g])
+(cf::CostFunction, g::OptTerm) = push!(cf.Terms,g)
