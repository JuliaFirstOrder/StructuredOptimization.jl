__precompile__()

module StructuredOptimization

using AbstractOperators
using AbstractOperators.BlockArrays
using ProximalOperators
using ProximalAlgorithms

import ProximalOperators:
			  is_affine,
			  is_cone,
			  is_convex,
			  is_generalized_quadratic,
			  is_prox_accurate,
			  is_quadratic,
			  is_separable,
			  is_set,
			  is_singleton,
			  is_smooth,
			  is_strongly_convex

include("calculus/precomposeNonlinear.jl")
include("syntax/syntax.jl")
include("solvers/solvers.jl")

export @minimize 

macro minimize(cf::Union{Expr, Symbol})
	cost = esc(cf)
	return :(solve(problem($(cost)), default_solver()))
end

macro minimize(cf::Union{Expr, Symbol}, s::Symbol, cstr::Union{Expr, Symbol})
	cost = esc(cf)
	if s == :(st)
		constraints = esc(cstr)
		return :(solve(problem($(cost), $(constraints)), default_solver()))
	elseif s == :(with)
		solver = esc(cstr)
		return :(solve(problem($(cost)), $(solver)))
	else
		error("wrong symbol after cost function! use `st` or `with`")
	end
end

macro minimize(cf::Union{Expr, Symbol}, s::Symbol, cstr::Union{Expr, Symbol}, w::Symbol, slv::Union{Expr, Symbol})
	cost = esc(cf)
	s != :(st) && error("wrong symbol after cost function! use `st`")
	constraints = esc(cstr)
	w != :(with) && error("wrong symbol after constraints! use `with`")
	solver = esc(slv)
	return :(solve(problem($(cost), $(constraints)), $(solver)))
end

end
