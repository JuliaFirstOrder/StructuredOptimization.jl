# check terms are proximable
# is_proximable{N}(xAll::NTuple{N,Variable}, t::Term) = is_AAc_diagonal(t)
#
# function is_proximable{N,M}(xAll::NTuple{N,Variable}, t::NTuple{M,Term})
# 	all(is_AAc_diagonal(t)) == false && return false
# 	var_appears = zeros(Int,N)
# 	is_there_a_GetIndex = false
# 	for i in eachindex(xAll)
# 		for ti in t
# 			if xAll[i] in variables(ti)
# 				# special case if operator is GetIndex,
# 				# there can be multiple terms with same varibles
# 				if typeof(operator(ti)) <: GetIndex
# 					if is_there_a_GetIndex == false var_appears[i] += 1 end
# 					is_there_a_GetIndex = true
# 				else
# 					var_appears[i] += 1
# 				end
# 			end
# 		end
# 	end
# 	return all(var_appears .<= 1) ? true : false
# end

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
		if length([t for t in terms if v in variables(t)]) != 1
			return false
		end
	end
	# NOTE: I see why GetIndex requires a special case. However, it is a more
	# general case than just GetIndex, and I would postpone its implementation,
	# unless we have a very concrete and important example where this is
	# strictly required...
	return true
end
