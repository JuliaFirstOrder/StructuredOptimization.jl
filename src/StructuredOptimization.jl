__precompile__()

module StructuredOptimization

using LinearAlgebra
using AbstractOperators
using ProximalOperators
using ProximalAlgorithms
using RecursiveArrayTools

include("syntax/syntax.jl")
include("calculus/precomposeNonlinear.jl") #TODO move to ProximalOperators?
include("solvers/solvers.jl")

end
