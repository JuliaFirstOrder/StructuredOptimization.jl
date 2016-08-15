module RegLS

include("prox.jl")
include("LBFGS.jl")
include("utils.jl")
include("ista.jl")
include("fista.jl")
include("zerofpr.jl")

solve = zerofpr

end
