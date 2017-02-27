#split cost funtion into smooth and non-smooth variables 

function split(cf::CostFunction)

	smooth = Array{OptTerm, 1}(0)
	nonsmooth = Array{OptTerm, 1}(0)
	#collect smooth terms (currently must be a ls)
	for i in eachindex(cf.Terms)
		if typeof(cf.Terms[i]) <: RegLS.SmoothTerm
			push!(smooth,cf.Terms[i])
		elseif typeof(cf.Terms[i]) <: RegLS.NonSmoothTerm
			push!(nonsmooth,cf.Terms[i])
		end
	end
	return smooth, nonsmooth

end
