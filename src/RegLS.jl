__precompile__()

module RegLS

using ProximalOperators
import ProximalOperators: RealOrComplex

import Base: in, *, <=, ==, sum, isempty, rank, norm
import ProximalOperators: is_convex

include("utils.jl")
include("operators.jl")
include("variables.jl")
include("affine.jl")
include("functions.jl")
include("composite.jl")
include("solvers.jl")
include("problems.jl")

end
