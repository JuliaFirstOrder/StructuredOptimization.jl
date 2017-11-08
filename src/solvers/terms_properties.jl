is_proximable(term::Term) = is_AAc_diagonal(term)

function is_proximable(terms::Tuple)
	# Check that each term is proximable
	if any(is_proximable.(terms) .== false)
		return false
	end
	# Construct the set of occurring variables
	vars = Set()
	for term in terms
		union!(vars, variables(term))
	end
	# Check that each variable occurs in only one term
	for v in vars
		tv = [t for t in terms if v in variables(t)]
		if length(tv) != 1
			#TODO make this more general
			if all( (<:).(typeof.(operator.(tv)), GetIndex) )
				return true
			else
				return false
			end
		end
	end
	# NOTE: I see why GetIndex requires a special case. However, it is a more
	# general case than just GetIndex, and I would postpone its implementation,
	# unless we have a very concrete and important example where this is
	# strictly required...
	# I agree... but we have this in the Audio Declipping demo!
	return true
end
