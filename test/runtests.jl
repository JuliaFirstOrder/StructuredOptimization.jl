using RegLS
using AbstractOperators
using ProximalOperators
using Base.Test
using Base.Profile

srand(0)

### Passing:

include("test_deep.jl")
include("test_functions.jl")
include("test_variables_expressions.jl")
include("test_terms.jl")
include("test_problem.jl")
include("test_solvers.jl")

### Yet to be checked:

# include("test_lbfgs.jl")
# include("test_lbfgs_larger.jl")
# include("test_costfunction.jl")
# include("test_matcomp.jl")
# include("test_svm.jl")
# include("test_minimize.jl")
# include("test_solvers_lasso.jl")
# include("test_solvers.jl")
