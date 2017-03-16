function split(cf::CostFunction)

	nonsmooth     = Array{NonSmoothTerm, 1}(0) # terms with non simple-operator
	quadratic     = Array{QuadraticTerm, 1}(0) # quadratic terms
	smooth        = Array{SmoothTerm,    1}(0) # quadratic terms
	nonsmoothprox = Array{NonSmoothTerm, 1}(0) # nonsmooth proximable terms

	for t in cf.Terms
		if isAbsorbable(t.A) && typeof(t)<:NonSmoothTerm
			push!(nonsmoothprox,t)
		else
			if typeof(t) <: NonSmoothTerm
				push!(nonsmooth,t)
			elseif typeof(t) <: SmoothTerm
				typeof(t) <: QuadraticTerm ? push!(quadratic,t) : push!(smooth,t)
			end
		end
	end

	return nonsmooth, smooth, quadratic, nonsmoothprox

end

