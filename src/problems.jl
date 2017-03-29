export problem, minimize

include("problems/split.jl")
include("problems/mergeProx.jl")
include("problems/mergeSmooth.jl")
include("problems/primal.jl")
#include("problems/dual.jl")

problem(h::CostFunction) = problem(h, Array{CostFunction,1}(0))
problem(h::CostFunction, cstr::CostFunction) = problem(h, [cstr])

function problem{T<:CostFunction}(cf::CostFunction, cstr::Array{T,1} )
	#add constraints to cost function
	for c in cstr
		cf += c
	end
	smooth, proximable, nonsmooth = split(cf)

	p = mergeProx(variable(cf), proximable)
	s = mergeSmooth(variable(cf), smooth)

	if isempty(nonsmooth)
		return Primal(s,p)
	else
		error("dual or smooth not implemented yet")
	end


end

minimize(h::CostFunction, args...) = minimize(h, Array{CostFunction,1}(0), args...)
minimize(h::CostFunction, cstr::CostFunction, args...) = minimize(h, [cstr], args...)

function minimize{T<:CostFunction}(cf::CostFunction, cstr::Array{T,1}, slv::Solver=ZeroFPR())
	P = problem(cf,cstr)
	return solve(P,slv)
end

		














