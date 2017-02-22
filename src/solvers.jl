abstract Solver
abstract ForwardBackwardSolver <: Solver

export
  PG,
  FPG,
  ZeroFPR

include("solvers/pg.jl")
include("solvers/zerofpr.jl")

include("solvers/utils.jl")

export solve, solve!

solve(A, b::AbstractArray, g::ProximableFunction) =
	solve(A, b::AbstractArray, g::ProximableFunction, zeros(eltype(b), size(A,2)))

solve(A::Union{AbstractArray,LinearOp}, b::AbstractArray, g::ProximableFunction, x::Array) =
  solve(A, b, g, x, ZeroFPR())

function solve(A::Union{AbstractArray,LinearOp}, b::AbstractArray, g::ProximableFunction, x0::AbstractArray, args...)
	x = deepcopy(x0) # copy initial conditions
	x, slv = solve!(A, b, g, x, args...)
	return x, slv
end
#
## when linear operator is composed with g (instead of the least squares term)
## then matrix/linear operator arguments is passed after g in the arguments list
## in this case the dual problem is solved, then the primal solution is recovered
#
function solve(b::AbstractArray, g::ProximableFunction, A, args...)
	y, slv = solve(A', b, Conjugate(g), args...)
	return -A'*y+b, slv, y
end
