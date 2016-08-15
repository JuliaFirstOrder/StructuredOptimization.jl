VERSION >= v"0.4.0-dev+6521" && __precompile__()

module RegLS

include("prox.jl")
include("utils.jl")
include("LBFGS.jl")
include("ista.jl")
include("fista.jl")
include("zerofpr.jl")

using Reexport 
@reexport using .Prox, .Ista, .Fista, .Zerofpr

solve = zerofpr
export solve
end
