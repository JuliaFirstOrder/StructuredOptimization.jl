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

#TODO move this somewhere else
function get_prox(x::Array{OptVar}, nonsmoothprox::Array{NonSmoothTerm,1})
	if length(x) == 1 #no block variables
		return get_prox(x[1],nonsmoothprox)
	else
		p = Array{ProximableFunction,1}(length(x))
		gxi   = Array{NonSmoothTerm,1}()
		for i in eachindex(x)
			for ns in nonsmoothprox
				if variable(ns) == x[i]
					push!(gxi, ns) 
				end
			end
			p[i] = get_prox(x[i],gxi)
			gxi = Array{NonSmoothTerm,1}()
		end
		return SeparableSum(p)
	end
end

function get_prox(x::OptVar, nonsmoothprox::Array{NonSmoothTerm,1} )
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
