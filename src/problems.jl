export minimize, minimize!

include("problems/split.jl")
include("problems/primal.jl")

#TODO implement also minimize!
minimize(h::OptTerm, args...) = minimize(CostFunction([h]),args...)

minimize{T<:OptTerm}(cf::CostFunction, cstr::Array{T}) = minimize(cf,cstr...)

function minimize(cf::CostFunction, cstr::Vararg{OptTerm})
	fi, smooth, nonsmooth = split(cf::CostFunction)
	for i = 1:length(cstr)
		push!(nonsmooth,cstr[i])
	end
	if length(fi) >  1 error("problem not supported") end
	if isempty(fi) 
		if isempty(smooth)
			error("only non smooth terms are present we should smooth a term")   
			# here we should have the call problem(nonsmooth)
		else
			push!(fi,pop!(smooth))
		end
	end 
	P = problem(fi[1], smooth, nonsmooth) # if fi[1] <: SmoothTerm creates Primal else Dual
	                                      # in Dual if  isempty(smooth) call problem(nonsmooth)
	solve(P,ZeroFPR())
end

		














