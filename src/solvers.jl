abstract Solver
abstract ForwardBackwardSolver <: Solver

export PG,
       FPG,
       ZeroFPR

include("solvers/pg.jl")
include("solvers/zerofpr.jl")

include("solvers/utils.jl")
include("solvers/lbfgs.jl")
