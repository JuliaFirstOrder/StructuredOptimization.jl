function mergeProx(x::Array{OptVar}, nonsmoothprox::Array{NonSmoothTerm,1})
	if length(x) == 1 #no block variables
		return mergeProx(x[1],nonsmoothprox)
	else
		p = Array{ProximableFunction,1}(length(x))
		gxi   = Array{NonSmoothTerm,1}()
		for i in eachindex(x)
			for ns in nonsmoothprox
				if variable(ns) == x[i]
					push!(gxi, ns) 
				end
			end
			p[i] = mergeProx(x[i],gxi)
			gxi = Array{NonSmoothTerm,1}()
		end
		return SeparableSum(p)
	end
end

function mergeProx(x::OptVar, nonsmoothprox::Array{NonSmoothTerm,1} )
	if length(nonsmoothprox) <= 1 
		if isempty(nonsmoothprox)
			p = IndFree()
		else
			p = absorbOp(nonsmoothprox[1].A, get_prox(nonsmoothprox[1]) )
		end
	else
		error("sliced separable sum not implemented, or there are multiple proximable terms with the same variables i.e. separable sum not possible!")
	end
end
