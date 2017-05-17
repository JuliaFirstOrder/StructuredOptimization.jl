export smooth

function smooth(cf::CompositeFunction, gamma::Real=1.)
	fs = Vector{ProximableFunction}(0)
	for f in terms(cf)
		s = MoreauEnvelope(gamma, f)
		push!(fs, s)
	end
	CompositeFunction(variable(cf), fs, affine(cf))
end
