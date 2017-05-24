using RegLS
using ProximalOperators
using Base.Test
using Base.Profile

srand(0)

### Passing:

include("test_operators.jl")
include("test_operators_calculus.jl")
include("test_functions.jl")
include("test_variables_expressions.jl")
include("test_terms.jl")
include("test_problem.jl")

### To be re-included:

# include("test_lbfgs.jl")
# include("test_lbfgs_larger.jl")

### Yet to be checked:

# include("test_costfunction.jl")
# include("test_matcomp.jl")
# include("test_svm.jl")
# include("test_minimize.jl")
# include("test_solvers_lasso.jl")

### This should be correct:

#include("test_solvers.jl")

### To be removed:

# include("test_linOp_old.jl")
