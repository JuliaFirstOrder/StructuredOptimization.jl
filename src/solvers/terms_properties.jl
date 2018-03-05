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
            if all( is_sliced.(tv) ) && all( is_proximable.(tv) )
				return true
			else
				return false
			end
		end
	end
	return true
end
