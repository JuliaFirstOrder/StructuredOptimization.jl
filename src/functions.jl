abstract OptTerm
abstract SmoothTerm    <: OptTerm
abstract NonSmoothTerm <: OptTerm

include("functions/absorb_merge.jl")
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

operator(h::OptTerm) = h.A
variable(h::OptTerm) = h.A.x

function Base.show(io::IO, f::OptTerm)
  println(io, "description : ", fun_name(f))
  println(io, "operator    : ", fun_name(f.A))
  println(io, "parameters  : ", fun_par(f))
end





