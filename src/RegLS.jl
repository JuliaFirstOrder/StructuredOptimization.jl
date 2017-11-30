__precompile__()

module RegLS

using AbstractOperators, ProximalOperators
using AbstractOperators: BlockArrays
import ProximalOperators: RealOrComplex,
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

#include("ProxStuff/functions/sqrHingeLoss.jl") #to MOVE to ProximalOperators
#include("ProxStuff/functions/crossEntropy.jl") #to MOVE to ProximalOperators
include("ProxStuff/calculus/precomposeNonlinear.jl") #to MOVE to ProximalOperators

include("syntax/syntax.jl")
include("solvers/solvers.jl")

export @minimize, st, with

immutable st end
immutable with end

macro minimize(cf::Union{Expr,Symbol})
	cost = esc(cf)
	return :(     solve!( problem($(cost)), default_slv())    )
end

macro minimize(cf::Union{Expr,Symbol}, s::Symbol, cstr::Union{Expr,Symbol})
	cost = esc(cf)
	if s == :(st)
		constraints = esc(cstr)
		return :(     solve!( problem($(cost),$(constraints)), default_slv())    )
	elseif s == :(with)
		solver = esc(cstr)
		return :(     solve!( problem($(cost)), $(solver))    )
	else
		error("wrong symbol after cost function! use `st` or `with`")
	end
end

macro minimize(cf::Union{Expr,Symbol}, s::Symbol, cstr::Union{Expr,Symbol}, w::Symbol, slv::Union{Symbol,Expr})
	cost = esc(cf)
	s != :(st) && error("wrong symbol after cost function! use `st`")
	constraints = esc(cstr)
	w != :(with) && error("wrong symbol after constraints! use `with`")
	solver = esc(slv)
	return :(     solve!( problem($(cost),$(constraints)), $(solver))    )
end

end
