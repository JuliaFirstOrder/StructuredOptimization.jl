function mergeProx(x_sorted::Vector{OptVar}, cf::CostFunction)
	if length(x_sorted) == 1 #no block variables
		return mergeProx(terms(cf), affine(cf))
	else
		p = Vector{ProximableFunction}(length(x_sorted))
		fi   = Vector{ExtendedRealValuedFunction}()
		Ai   = Vector{AffineOperator}()
		for i in eachindex(x_sorted)
			for ii in eachindex(affine(cf))
				if variable(affine(cf)[ii])[1] == x_sorted[i]
					push!(fi,  terms(cf)[ii]) 
					push!(Ai, affine(cf)[ii]) 
				end
			end
			p[i] = mergeProx(fi,Ai)
		
			fi   = Vector{ExtendedRealValuedFunction}() #reinitialize the arrays
			Ai   = Vector{AffineOperator}()
		end
		return SeparableSum(p)
	end
end

function mergeProx(f::Vector{ExtendedRealValuedFunction}, affOps::Vector{AffineOperator})
	if length(f) <= 1 
		if isempty(f)
			p = IndFree()
		else
			p = absorbOp(affOps[1], get_prox(f[1]) )
		end
	else
		error("sliced separable sum not implemented, or there are multiple proximable terms with the same variables i.e. separable sum not possible!")
	end
end
