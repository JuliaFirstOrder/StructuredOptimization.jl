abstract OptTerm
abstract SmoothTerm    <: OptTerm
abstract NonSmoothTerm <: OptTerm

include("functions/get_prox.jl")
include("functions/LeastSquares.jl")
include("functions/HingeLoss.jl")
include("functions/Norm.jl")
include("functions/Box.jl")

type CostFunction
	Terms::Array{OptTerm, 1}
end

+(h::OptTerm, g::OptTerm) = CostFunction([h,g])

function +(cf::CostFunction, g::OptTerm) 
	push!(cf.Terms,g)
	return cf
end

#returns the number of blocks of variables
function blkLength(h::OptTerm)
	blkLength(h.A.x)
end

operator(h::OptTerm) = h.A
variable(h::OptTerm) = h.A.x






