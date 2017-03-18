export minimize #, minimize!

include("problems/split.jl")
include("problems/mergeProx.jl")
include("problems/mergeSmooth.jl")
include("problems/primal.jl")
include("problems/dual.jl")

#TODO implement also minimize!
minimize(h::OptTerm, args...)                         = minimize(CostFunction([variable(h)] ,[h]),args...)
minimize(h::CostFunction, args...)                    = minimize(h, Array{OptTerm,1}(0), args...)
minimize(h::CostFunction, cstr::OptTerm, args...)     = minimize(h, [cstr], args...)

function minimize{T<:OptTerm}(cf::CostFunction, cstr::Array{T,1}, slv::Solver = ZeroFPR() )
	#add constraints to cost function
	for c in cstr
		cf += c
	end
	nonsmooth, smooth, quadratic, nonsmoothprox = split(cf)
	# here these are all arrays which should be merged into a single term
	prm = sortperm(variable(cf),by = object_id)
	x_sorted = sort(variable(cf),by = object_id)
	# for example for the nonsmoothprox these are merge in such a way:
	prox = mergeProx(x_sorted,nonsmoothprox) #where prox is a ProximableFunction

	#merging quadratic is still not done
	quad = mergeSmooth(x_sorted, quadratic)

	x, slv = solve(quad, prox, slv) #this should change with the new call to solve

	if prm == [1]
		return x, slv
	else
		x .=x[prm]
		return x, slv
	end

end

		














