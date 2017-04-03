function split(cf::CostFunction)

	nA, nf = Vector{AffineOperator}(0), Vector{ExtendedRealValuedFunction}(0) 
	# non-smooth terms with non simple-operator
	sA, sf = Vector{AffineOperator}(0), Vector{ExtendedRealValuedFunction}(0) 
	# smooth terms 
	pA, pf = Vector{AffineOperator}(0), Vector{ExtendedRealValuedFunction}(0) 
	# nonsmooth proximable terms

	for i in eachindex(cf.A)
		if isAbsorbable(operator(cf.A[i])) && typeof(cf.f[i])<:NonSmoothFunction
			push!(pf,cf.f[i])
			push!(pA,cf.A[i])
		else
			if     typeof(cf.f[i]) <: NonSmoothFunction
				push!(nf,cf.f[i])
				push!(nA,cf.A[i])
			elseif typeof(cf.f[i]) <: SmoothFunction
				push!(sf,cf.f[i])
				push!(sA,cf.A[i])
			end
		end
	end

	proximable =  CostFunction(variable(cf), pf, pA )
	smooth     =  CostFunction(variable(cf), sf, sA )
	nonsmooth  =  CostFunction(variable(cf), nf, nA )

	return smooth, proximable, nonsmooth 

end

split_quadratic(cf::CostFunction) = 
CostFunction(variable(cf), 
	      terms(cf)[isquadratic.(terms(cf))], 
	     affine(cf)[isquadratic.(terms(cf))]),
CostFunction(variable(cf), 
	      terms(cf)[!isquadratic.(terms(cf))], 
	     affine(cf)[!isquadratic.(terms(cf))])


