include("problems/extract.jl")
include("problems/terms_properties.jl")
include("problems/split.jl")
include("problems/solve.jl")

export problem, @minimize, st, with

function problem(terms::Vararg)
	cf = ()
	for i = 1:length(terms)
		cf = (cf...,terms[i]...)
	end
	return cf
end

immutable st end
immutable with end

macro minimize(cf::Expr) 
	cost = esc(cf)
	return :(     solve!( problem($(cost)), default_slv())    )
end

macro minimize(cf::Expr, s::Symbol, cstr::Expr) 
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

macro minimize(cf::Expr, s::Symbol, cstr::Expr, w::Symbol, slv::Expr) 
	cost = esc(cf)
	s != :(st) && error("wrong symbol after cost function! use `st`") 
	constraints = esc(cstr)
	w != :(with) && error("wrong symbol after constraints! use `with`") 
	solver = esc(slv)
	return :(     solve!( problem($(cost),$(constraints)), $(solver))    )
end

