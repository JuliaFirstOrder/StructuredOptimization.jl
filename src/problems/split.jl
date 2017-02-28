#split cost funtion into smooth and non-smooth variables 

function split(cf::CostFunction)

	fj = Array{OptTerm, 1}(0)
	smooth = Array{OptTerm, 1}(0)
	nonsmooth = Array{OptTerm, 1}(0)
	#collect smooth terms (currently must be a ls)
	for t in cf.Terms
		if isMergeable(t.A)
			if     typeof(t) <: LeastSquares && isEye(t.A)
				push!(smooth,t)
			elseif typeof(t) <: NonSmoothTerm
				push!(nonsmooth,t)
			end
		else
			push!(fj,t) 
		end
	end
	return fj, smooth, nonsmooth

end
