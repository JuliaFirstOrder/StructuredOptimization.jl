#split cost funtion into smooth and non-smooth variables 

function split(cf::CostFunction)

	g         = Array{OptTerm, 1}(0)
	smooth    = Array{OptTerm, 1}(0)
	nonsmooth = Array{OptTerm, 1}(0)

	for t in cf.Terms
		if isAbsorbable(t.A)
			if     typeof(t) <: LeastSquares && isDiagonal(t.A)
				push!(smooth,t)
			elseif typeof(t) <: NonSmoothTerm
				push!(nonsmooth,t)
			end
		else
			push!(g,t) 
		end
	end
	return g, smooth, nonsmooth

end
