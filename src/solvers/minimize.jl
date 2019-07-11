export @minimize

"""
    @minimize cost [st ctr] [with slv_opt]

Minimize a given problem with cost function `cost`, constraints `ctr` and solver options `slv_opt`.

# Example

```julia
julia> using StructuredOptimization

julia> A, b, x = randn(10,4), randn(10), Variable(4);

julia> @minimize ls(A*x-b) + 0.5*norm(x);

julia> ~x  # access array with solution

julia> @minimize ls(A*x-b) st x >= 0.;

julia> ~x  # access array with solution

julia> @minimize ls(A*x-b) st norm(x) == 2.0 with ForwardBackward(fast=true);

julia> ~x  # access array with solution
```

Returns as output a tuple containing the optimization variables and the number
of iterations spent by the solver algorithm.
"""
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
