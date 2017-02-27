export minimize, minimize!

include("problems/split.jl")
include("problems/primal.jl")

#TODO implement also minimize!
minimize(h::OptTerm, args...) = minimize(CostFunction([h]),args...)

minimize{T<:OptTerm}(cf::CostFunction, cstr::Array{T}) = minimize(cf,cstr...)

function minimize(cf::CostFunction, cstr::Vararg{OptTerm})
	smooth,nonsmooth = split(cf::CostFunction)
	for i = 1:length(cstr)
		push!(nonsmooth,cstr[i])
	end
	A,g = Primal(smooth,nonsmooth)
	solve(A,g,ZeroFPR())
end

		














