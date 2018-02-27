export @minimize 

"""
`@minimize cost [st ctr] [with slv_opt]`

Minimize a given problem with cost function `cost`, constraints `ctr` and solver options `slv_opt`. 

# Example

```julia
julia> A, b = randn(10,4), randn(10);

julia> @minimize ls(A*x-b) + 0.5*norm(x);
    it |      gamma |        fpr |        tau |        FBE |
 ------|------------|------------|------------|------------|
     1 | 2.9152e-02 | 2.7656e+00 | 1.0000e+00 | 5.5181e+00 |
     9 | 2.9152e-02 | 9.9682e-05 | 1.0000e+00 | 4.4086e+00 |

julia> @minimize ls(A*x-b) st x >= 0.;
    it |      gamma |        fpr |        tau |        FBE |
 ------|------------|------------|------------|------------|
     1 | 5.8304e-02 | 1.0068e+00 | 1.0000e+00 | 6.6282e+00 |
     3 | 5.8304e-02 | 9.5210e-16 | 1.0000e+00 | 6.5654e+00 |

julia> it, slv = @minimize ls(A*x-b) st norm(x) == 2.0 with PG(maxit = 5);
    it |      gamma |        fpr |
 ------|------------|------------|
     1 | 6.1373e-02 | 2.2090e+01 |
     5 | 3.0686e-02 | 5.5190e-01 |

```

Returns as output a tuple containing the number of iterations and the constructed solver.

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
