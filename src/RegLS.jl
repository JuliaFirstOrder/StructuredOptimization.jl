__precompile__()

module RegLS

using ProximalOperators
import ProximalOperators: RealOrComplex

# import Base: in, +, *, <=, ==, sum, isempty, rank, norm

include("utils.jl")
include("operators.jl")
include("functions.jl")
include("variables.jl")
include("affine.jl")
include("term.jl")
include("solvers.jl")
include("problem.jl")

end
