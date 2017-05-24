# export smooth
#
# function smooth(cf::Term, gamma::Real=1.)
# 	fs = Vector{ProximableFunction}(0)
# 	for f in terms(cf)
# 		s = MoreauEnvelope(gamma, f)
# 		push!(fs, s)
# 	end
# 	Term(variable(cf), fs, affine(cf))
# end
