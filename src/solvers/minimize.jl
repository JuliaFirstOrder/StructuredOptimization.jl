export @minimize

"""
    @minimize cost [st ctr] [with slv_opt]

Minimize a given problem with cost function `cost`, constraints `ctr` and solver options `slv_opt`.

# Example

```julia
julia> using StructuredOptimization, ProximalAlgorithms

julia> A, b, x = randn(10,4), randn(10), Variable(4);

julia> @minimize ls(A*x-b) + 0.5*norm(x);

julia> ~x
4-element Array{Float64,1}:
 -0.3606117955906451
  0.21231949195945288
  0.022423882830584645
  0.19560804140788765

julia> @minimize ls(A*x-b) st x >= 0.;

julia> ~x
4-element Array{Float64,1}:
 0.0
 0.2682589660003247
 0.0
 0.2297141914608737

julia> @minimize ls(A*x-b) st norm(x) == 2.0 with ProximalAlgorithms.ForwardBackward(fast=true);

julia> ~x
4-element Array{Float64,1}:
 -0.7217037749126615
  0.9835726063727642
  0.634450485411672
  1.452308910263514

julia>
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
