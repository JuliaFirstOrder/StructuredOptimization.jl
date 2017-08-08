
# check terms are proximable
is_proximable{N}(xAll::NTuple{N,Variable}, t::Term) = is_AAc_diagonal(t)

function is_proximable{N,M}(xAll::NTuple{N,Variable}, t::NTuple{M,Term})
	all(is_AAc_diagonal(t)) == false && return false
	var_appears = zeros(Int,N)
	is_there_a_GetIndex = false
	for i in eachindex(xAll)
		for ti in t 
			if xAll[i] in variables(ti) 
				# special case if operator is GetIndex, 
				# there can be multiple terms with same varibles
				if typeof(operator(ti)) <: GetIndex 
					if is_there_a_GetIndex == false var_appears[i] += 1 end 
					is_there_a_GetIndex = true
				else
					var_appears[i] += 1
				end
			end
		end
	end
	return all(var_appears .<= 1) ? true : false
end

