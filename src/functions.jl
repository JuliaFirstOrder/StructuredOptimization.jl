abstract OptTerm
abstract SmoothTerm      <: OptTerm
abstract QuadraticTerm   <: SmoothTerm
abstract NonSmoothTerm   <: OptTerm

include("functions/CostFunction.jl")
include("functions/absorb_merge.jl")
include("functions/LeastSquares.jl")
include("functions/HingeLoss.jl")
include("functions/Norm.jl")
include("functions/Box.jl")

operator(h::OptTerm) = h.A
variable(h::OptTerm) = variable(h.A)
variable(cf::CostFunction) = cf.x

function Base.show(io::IO, f::OptTerm)
  println(io, "description : ", fun_name(f))
  println(io, "operator    : ", fun_name(f.A))
  println(io, "parameters  : ", fun_par(f))
end





