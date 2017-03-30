abstract ExtendedRealValuedFunction
abstract SmoothFunction      <: ExtendedRealValuedFunction
abstract NonSmoothFunction   <: ExtendedRealValuedFunction

include("functions/CostFunction.jl")
include("functions/absorb_merge.jl")
include("functions/LeastSquares.jl")
include("functions/Norm.jl")
include("functions/Box.jl")
include("functions/HingeLoss.jl")

gradient(f::ExtendedRealValuedFunction) = error("gradient not implemented for $f")






