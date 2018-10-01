__precompile__()

module StructuredOptimization

using LinearAlgebra
using AbstractOperators
using AbstractOperators.BlockArrays
using ProximalOperators
using ProximalAlgorithms

include("syntax/syntax.jl")
include("calculus/precomposeNonlinear.jl") #TODO move to ProximalOperators?
include("solvers/solvers.jl")

end
