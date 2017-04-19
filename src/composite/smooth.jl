export smooth

function smooth(cf::CompositeFunction, gamma0::Real=1.)
	f = Vector{ExtendedRealValuedFunction}(0)
	for fi in terms(cf)
		if typeof(fi)<:SmoothFunction
			push!(f,fi)
		elseif typeof(fi)<:NonSmoothFunction
			p = get_prox(fi)
			fs = MoreauEnvelope(gamma0*lambda(fi),p) 
			push!(f,fs)
		end
	end
	CompositeFunction(variable(cf),f,affine(cf))
end

