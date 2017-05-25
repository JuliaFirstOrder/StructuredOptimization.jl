__precompile__()

module RegLS

using ProximalOperators
import ProximalOperators: RealOrComplex

include("deep.jl")
include("operators.jl")
include("functions.jl")
include("syntax.jl")
include("solvers.jl")
include("problem.jl")

end
