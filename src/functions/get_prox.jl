
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
				mergeProx(p, reg[1].lambda, reg[1].A)  
			else
				error("cannot merge $(typeof(reg[1])) with $(typeof(nonsmooth[1].A))")
			end
		end
	else
		error("sliced separable sum not implemented yet")
	end
end

# absorb linear operator into proximable operator
absorbOp(A::LinearOperator, p::ProximableFunction) = absorbOp(  A, p,  0)
absorbOp(A::Affine  , p::ProximableFunction) = absorbOp(A.A, p,A.b)

absorbOp{L <:IdentityOperator}(A::L, p::ProximableFunction, b) = b == 0 ? p : Precompose(p, 1, b)
absorbOp{L <:DiagonalOperator}(A::L, p::ProximableFunction, b) = Precompose(p, A.d, b)

# merge Proximal operators
mergeProx{T<:AffineOperator}(p::ProximableFunction, lambda, A::T) =
Regularize(p, lambda, -A.b)

mergeProx{T<:LinearOperator}(p::ProximableFunction, lambda, A::T) =
Regularize(p, lambda, 0.)





