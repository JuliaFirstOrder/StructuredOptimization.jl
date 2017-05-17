export problem, minimize

include("problems/primal.jl")
# include("problems/dual.jl")
# include("problems/smoothed.jl")

problem(h::CompositeFunction) = problem(h, Array{CompositeFunction,1}(0))
problem(h::CompositeFunction, cstr::CompositeFunction) = problem(h, [cstr])

function problem{T<:CompositeFunction}(cf::CompositeFunction, cstr::Array{T,1} )
	#add constraints to cost function
	for c in cstr
		cf += c
	end
	smooth, proximable, nonsmooth = split(cf)

	proximable = mergeProx(variable(cf), proximable)
	smooth     = sort_and_expand(variable(cf), smooth    )
	nonsmooth  = sort_and_expand(variable(cf), nonsmooth    )

	if isempty(nonsmooth)
		return Primal(smooth,proximable)
	else
		if isempty(proximable) && length(terms(nonsmooth)) == 1
			return Dual(smooth,nonsmooth)
		else
			return SmoothedPrimal(smooth,nonsmooth,proximable)
		end
	end
end

minimize(h::CompositeFunction, args...) = minimize(h, Array{CompositeFunction, 1}(0), args...)
minimize(h::CompositeFunction, cstr::CompositeFunction, args...) = minimize(h, [cstr], args...)

function minimize{T <: CompositeFunction}(cf::CompositeFunction, cstr::Array{T,1}, slv::Solver=ZeroFPR())
	P = problem(cf, cstr)
	return solve(P, slv)
end
