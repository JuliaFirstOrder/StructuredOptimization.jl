abstract Solver
abstract ForwardBackwardSolver <: Solver

export PG,
       FPG,
       ZeroFPR

include("solvers/utils.jl")
include("solvers/lbfgs.jl")

include("solvers/pg.jl")
include("solvers/fpg.jl")
include("solvers/zerofpr.jl")

return
