
get_prox(A::LinearOp, g::ProximableFunction) = get_prox(A, g, 0)

get_prox{L <:IdentityOperator}(A::L, g::ProximableFunction, b) = b == 0 ? g : Precompose(g, 1, b)
get_prox{L <:DiagonalOperator}(A::L, g::ProximableFunction, b) = Precompose(g, A.d, b)

get_prox(A::Affine, g::ProximableFunction) = get_prox(A.A,g,A.b)

function get_prox(x::OptVar, nonsmooth::Array{OptTerm,1}, reg::Array{OptTerm,1} )
	if length(nonsmooth) <= 1  
		if isempty(nonsmooth)
			g = IndFree()
		else
			if x == nonsmooth[1].A.x
				if ismergable(nonsmooth[1].A)
					g = get_prox(nonsmooth[1].A, get_prox(nonsmooth[1]) )
				else
				error("cannot merge linear operator $(typeof(nonsmooth[1].A)) with proximal operator, try dual?")
				end
			else
				error("terms have different variables!")
			end
		end
	else
		error("sliced separable sum not implemented yet")
	end
	if isempty(reg)
		return g
	else
		length(reg) == 1 ? Regularize(g, reg[1].lambda, 0.) : error("")
	end
end

function get_prox(x::Array{OptVar}, nonsmooth::Array{OptTerm,1}, reg::Array{OptTerm,1})
	g = Array{ProximableFunction,1}(length(x))
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
		g[i] = get_prox(x[i],gxi,regxi)
		gxi = Array{OptTerm,1}()
		regxi = Array{OptTerm,1}()
	end
	return SeparableSum(g)

end

ismergable(A::Affine) = ismergable(A.A)
ismergable(A::LinearOp) = typeof(A) <: DiagonalOperator

