using RegLS
using ProximalOperators
using Base.Test
using Base.Profile

srand(0)

### Passing:

include("test_linOp.jl")
include("test_linOp_calculus.jl")
include("test_lbfgs.jl")
include("test_lbfgs_larger.jl")
include("test_functions.jl")

### Yet to be checked:

# include("test_affine_and_variables.jl")
# include("test_costfunction.jl")
# include("test_matcomp.jl")
# include("test_svm.jl")
# include("test_minimize.jl")
# include("test_solvers_lasso.jl")

### This should be correct:

#include("test_solvers.jl")
