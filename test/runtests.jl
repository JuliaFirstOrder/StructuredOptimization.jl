using RegLS
using ProximalOperators
using Base.Test
using Base.Profile
using MathProgBase
using Ipopt

srand(1)

#include("test_linOp.jl")
include("test_functions.jl")
#include("test_minimize.jl")
#include("test_lbfgs.jl")
#include("test_lbfgs_larger.jl")
#include("test_solvers_lasso.jl")
#include("test_matcomp.jl")
#include("test_svm.jl")
