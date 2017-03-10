immutable Primal{A<:AffineOperator}
	Op::A
	p::ProximableFunction
end

function problem{T<:SmoothTerm}(g::T, smooth::Array{OptTerm,1}, nonsmooth::Array{OptTerm,1})

	x = variable(g) #extract variables
	p = get_prox(x, nonsmooth, smooth)
	 
	Primal(g.A, p)
end

solve(P::Primal, args...) = solve(P.Op, P.p, args...)

function get_prox(x::Array{OptVar}, nonsmooth::Array{OptTerm,1}, reg::Array{OptTerm,1})
	if length(x) == 1
		get_prox(x[1],nonsmooth,reg) #this happens with VCAT
	else
		p = Array{ProximableFunction,1}(length(x))
		gxi   = Array{OptTerm,1}()
		regxi = Array{OptTerm,1}()
		for i in eachindex(x)
			for ns in nonsmooth
				if ns.A.x == x[i]
					push!(gxi, ns) 
				end
			end
			for r in reg
				if r.A.x == x[i]
					push!(regxi, r) 
				end
			end
			p[i] = get_prox(x[i],gxi,regxi)
			gxi = Array{OptTerm,1}()
			regxi = Array{OptTerm,1}()
		end
		return SeparableSum(p)
	end
end

function get_prox(x::OptVar, nonsmooth::Array{OptTerm,1}, reg::Array{OptTerm,1} )
	if length(nonsmooth) <= 1 && length(reg) <= 1 
		if isempty(nonsmooth)
			p = IndFree()
		else
			if isAbsorbable(nonsmooth[1].A)
				p = absorbOp(nonsmooth[1].A, get_prox(nonsmooth[1]) )
			else
				error("cannot absorb linear operator $(typeof(nonsmooth[1].A))")
			end
		end
		if isempty(reg)
			return p
		else
			if reg[1].A == nonsmooth[1].A
				return  mergeProx(p, reg[1].lambda, reg[1].A)  
			else
				error("cannot merge $(typeof(reg[1])) with $(typeof(nonsmooth[1].A))")
			end
		end
	else
		error("sliced separable sum not implemented yet")
	end
end
