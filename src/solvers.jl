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

solve(A::AbstractArray, b::AbstractArray, g::ProximableFunction) =
solve(A, b, g, zeros(eltype(b), size(A,2)))

solve(A::AbstractArray, b::AbstractArray, g::ProximableFunction, x0::AbstractArray) = 
solve(A, b, g, x0, ZeroFPR() )

function solve(A::AbstractArray, b::AbstractArray, g::ProximableFunction, x0::AbstractArray, args...) 
	x = OptVar(deepcopy(x0))
	L = A*x+b
	solve!(L, g, args...)
end

#
## when linear operator is composed with g (instead of the least squares term)
## then matrix/linear operator arguments is passed after g in the arguments list
## in this case the dual problem is solved, then the primal solution is recovered
#
function solve{T <: AffineOperator}(A::T, args...)
	x0 = optArray(A)
	x, slv = solve!(A, args...)
	optArray!(A,x0)
	return x, slv
end

function solve(b::AbstractArray, g::ProximableFunction, A::AbstractArray, args...)
	y, slv = solve(A', b, Conjugate(g), args...)
	return -A'*y+b, slv, y
end

function solve!(g::ProximableFunction, A::AffineOperator, args...)
	At = A' + A.b
	y, slv = solve!(At, Conjugate(g), args...)
	return -At'*y, slv, y
end

function solve!(g::ProximableFunction, A::LinearOp, args...)
	y, slv = solve!(A', Conjugate(g), args...)
	return -A'*y, slv, y
end
