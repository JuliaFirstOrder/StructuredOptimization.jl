#split cost funtion into smooth and non-smooth variables 

function split(cf::CostFunction)

	g         = Array{OptTerm, 1}(0) # terms with non-simple linear operators
	smooth    = Array{OptTerm, 1}(0) # smooth terms
	nonsmooth = Array{OptTerm, 1}(0) # non smooth terms

	for t in cf.Terms
		if isAbsorbable(t.A)
			if     typeof(t) <: LinearLeastSquares && isDiagonal(t.A)
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
