__precompile__()

module RegLS

using AbstractOperators, ProximalOperators
import ProximalOperators: RealOrComplex

include("deep.jl")
include("functions.jl")
include("syntax.jl")
include("solvers.jl")
include("problem.jl")

end
