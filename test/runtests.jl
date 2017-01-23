using RegLS
using ProximalOperators
using Base.Test
using Base.Profile
using MathProgBase
using Ipopt

srand(0)

include("test_lbfgs.jl")
include("test_lbfgs_larger.jl")
include("test_lasso_random.jl")
include("test_lasso_complex.jl")
include("test_lasso_multidim.jl")
include("test_matcomp.jl")
include("test_multi_var.jl")
include("test_svm.jl")
