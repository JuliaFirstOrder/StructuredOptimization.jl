export problem, minimize

# include("problems/primal.jl")
# include("problems/dual.jl")
# include("problems/smoothed.jl")

# problem(h::Term) = problem(h, Array{Term,1}(0))
# problem(h::Term, cstr::Term) = problem(h, [cstr])

# function problem{T<:Term}(cf::Term, cstr::Array{T,1} )
# 	#add constraints to cost function
# 	for c in cstr
# 		cf += c
# 	end
# 	smooth, proximable, nonsmooth = split(cf)
#
# 	proximable = mergeProx(variable(cf), proximable)
# 	smooth     = sort_and_expand(variable(cf), smooth    )
# 	nonsmooth  = sort_and_expand(variable(cf), nonsmooth    )
#
# 	if isempty(nonsmooth)
# 		return Primal(smooth,proximable)
# 	else
# 		if isempty(proximable) && length(terms(nonsmooth)) == 1
# 			return Dual(smooth,nonsmooth)
# 		else
# 			return SmoothedPrimal(smooth,nonsmooth,proximable)
# 		end
# 	end
# end

function problem(terms...)
	# Concatenate all terms
	cf = vcat(terms...)
	# Build up the set of all occurring variables
	# V = Set()
	# for term in cf
	# 	union!(V, variables(term.A))
	# end
	# # Expand affine expressions in each term
	# for term in cf
	# 	for v in V
	# 		# TODO: do something
	# 	end
	# end
end

minimize(cf::Vararg{Term}; solver::Solver=default_solver()) =
	solve(problem(cf...), solver)

# minimize(h::Term, args...) = minimize(h, Array{Term, 1}(0), args...)
# minimize(h::Term, cstr::Term, args...) = minimize(h, [cstr], args...)

# function minimize{T <: Term}(cf::Term, cstr::Array{T,1}, slv::Solver=ZeroFPR())
	# P = problem(cf, cstr)
	# return solve(P, slv)
# end
