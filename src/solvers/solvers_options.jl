abstract type Solver end

abstract type ForwardBackwardSolver <: Solver end

export Solver, ForwardBackwardSolver

################################################################################
export PG, FPG

"""
`PG(;kwargs...)`

Creates an object `PG` containing the options of the Proximal Gradient solvers:

  * `gamma`, stepsize (default: unspecified, determined automatically)
  * `maxit`, maximum number of iteration (default: `10000`)
  * `tol`, halting tolerance on the fixed-point residual (default: `1e-4`)
  * `adaptive`, adaptively adjust `gamma` (default: `false` if `gamma` is provided)
  * `fast`, enables accelerated method (default: `false`)
  * `verbose`, verbosity level (default: `1`)
  * `verbose_freq`, verbosity frequency for `verbose = 1` (default: `100`)

"""
struct PG <: ForwardBackwardSolver
    kwargs::Iterators.Pairs
    function PG(; kwargs...)
        new(kwargs)
    end
end

"""
`FPG(;kwargs...)`

Same as `PG`, creates the options of the Fast Proximal Gradient solver. 

"""
function FPG(; kwargs...)
    return PG(; kwargs..., fast=true)
end

function build_iterator(x, solver::PG; kwargs...)
    x, ProximalAlgorithms.FBSIterator(~x; solver.kwargs..., kwargs...)
end

################################################################################
export ZeroFPR

"""
`ZeroFPR(;kwargs...)`

Creates an object `ZeroFPR` containing the options of the ZeroFPR solver:

  * `gamma`, stepsize (default: unspecified, determined automatically)
  * `maxit`, maximum number of iteration (default: `10000`)
  * `tol`, halting tolerance on the fixed-point residual (default: `1e-4`)
  * `adaptive`, adaptively adjust `gamma` (default: `false` if `gamma` is provided)
  * `fast`, enables accelerated method (default: `false`)
  * `verbose`, verbosity level (default: `1`)
  * `verbose_freq`, verbosity frequency for `verbose = 1` (default: `100`)
  * `memory`, memory of the `LBFGS` operator (default: `10` )

"""
struct ZeroFPR <: ForwardBackwardSolver
    kwargs::Iterators.Pairs
    function ZeroFPR(; kwargs...)
        new(kwargs)
    end
end

function build_iterator(x, solver::ZeroFPR; kwargs...)
    x, ProximalAlgorithms.ZeroFPRIterator(~x; solver.kwargs..., kwargs...)
end

################################################################################
export PANOC

"""
`ZeroFPR(;kwargs...)`

Creates an object `PANOC` containing the options of the PANOC solver:

  * `gamma`, stepsize (default: unspecified, determined automatically)
  * `maxit`, maximum number of iteration (default: `10000`)
  * `tol`, halting tolerance on the fixed-point residual (default: `1e-4`)
  * `adaptive`, adaptively adjust `gamma` (default: `false` if `gamma` is provided)
  * `fast`, enables accelerated method (default: `false`)
  * `verbose`, verbosity level (default: `1`)
  * `verbose_freq`, verbosity frequency for `verbose = 1` (default: `100`)
  * `memory`, memory of the `LBFGS` operator (default: `10` )

"""
struct PANOC <: ForwardBackwardSolver
    kwargs::Iterators.Pairs
    function PANOC(; kwargs...)
        new(kwargs)
    end
end

function build_iterator(x, solver::PANOC; kwargs...)
    x, ProximalAlgorithms.PANOCIterator(~x; solver.kwargs..., kwargs...)
end

default_solver = PANOC

################################################################################
