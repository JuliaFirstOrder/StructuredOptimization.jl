
function Primal(smooth::Array{OptTerm,1}, nonsmooth::Array{OptTerm,1})
	#analyze smooth terms
	if length(smooth) == 0 error("no smooth term is present!") end
	f, reg = get_mergeable(smooth)

	x = f.A.x #extract variables
	g = get_prox(x,nonsmooth,reg)
	 
	return f.A, g
	
end

function get_mergeable(smooth::Array{OptTerm,1})
	#choose the term with most block variables
	idx = findmax(blkLength.(smooth))
	if length(idx) == 1     #there is only one term with multiple var
		f = smooth[idx]
	else                
		idx = 0
		for i in eachindex(smooth) 
			idx += 1
			if isEye(smooth[i].A) == false 
				#get the term without identity operator
				#else get the last
				break
			end
		end
		f = smooth[idx]
	end
	reg = get_reg(smooth,idx)
	return f, reg
end

function get_reg(smooth::Array{OptTerm,1}, idx::Int64)
	reg = Array{OptTerm,1}()
	#other smooth terms must have identity and may be merged
	for i in eachindex(smooth) 
		if i != idx
			(typeof(smooth[i].A) <: IdentityOperator && 
                         typeof(smooth[i])   <: LeastSquares ) ? 
			push!(reg,smooth[i]) : error("cannot merge smooth term")
		end
	end
	return reg
end

