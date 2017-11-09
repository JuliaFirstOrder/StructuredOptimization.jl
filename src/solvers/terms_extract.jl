# returns all variables of a cost function, in terms of appearance
extract_variables(t::Term) = variables(t) 

function extract_variables{N}(t::NTuple{N,Term})
	x = variables.(t)
	xAll = x[1]
	for i = 2:length(x)
		for xi in x[i]
			if (xi in xAll) == false
				xAll = (xAll...,xi)
			end
		end
	end
	return xAll
end

# extract functions from terms
function extract_functions(t::Term)
	f = displacement(t.A) == 0 ? t.f : PrecomposeDiagonal(t.f, 1.0, displacement(t.A)) #for now I keep this
	f = t.lambda == 1. ? f : Postcompose(f, t.lambda)                                  #for now I keep this
	#TODO change this
	return f
end
extract_functions{N}(t::NTuple{N,Term}) = SeparableSum(extract_functions.(t))
extract_functions(t::Tuple{Term}) = extract_functions(t[1])

# extract operators from terms

# returns all operators with an order dictated by xAll

#single term, single variable
extract_operators(xAll::Tuple{Variable}, t::Term)  = operator(t)

extract_operators{N}(xAll::NTuple{N,Variable}, t::Term) = extract_operators(xAll, (t,))

#multiple terms, multiple variables
function extract_operators{N,M}(xAll::NTuple{N,Variable}, t::NTuple{M,Term})
	ops = ()
	for ti in t
		tex = expand(xAll,ti)
		ops = (ops...,sort_and_extract(xAll,tex))
	end
	return vcat(ops...)
end

function expand{N,T1,T2,T3}(xAll::NTuple{N,Variable}, t::Term{T1,T2,T3})
	xt   = variables(t)
	C    = codomainType(operator(t))
	size_out = size(operator(t),1)
	ex = t.A

	for x in xAll
		if !( x in variables(t) ) 
			ex += Zeros(eltype(~x),size(x),C,size_out)*x
		end
	end
	return Term{T1,T2,typeof(ex)}(t.lambda, t.f, ex)
end

sort_and_extract(xAll::Tuple{Variable}, t::Term) = operator(t)

function sort_and_extract{N}(xAll::NTuple{N,Variable}, t::Term)
	p = zeros(Int,N)
	xL = variables(t)
	for i in eachindex(xAll)
		p[i] = findfirst( xi -> xi == xAll[i], xL)
	end
	return operator(t)[p]
end

function extract_proximable(xAll::NTuple{N,Variable}, t::NTuple{M,Term}) where {N,M}
	fs = ()
	for x in xAll
		tx = () #terms containing x
		for ti in t
			if x in variables(ti)
				tx = (tx...,ti) #collect terms containing x
			end
		end
		if isempty(tx)
			fx = IndFree()
		elseif length(tx) == 1          #only one term per variable
			fx = extract_proximable(x,tx[1])
		else                            #multiple terms per variable
			                        #currently this happens only with GetIndex
		
			fxi,idxs = (),()
			for ti in tx
				fxi  = (fxi..., extract_functions(ti))
				idxs = (idxs...,operator(ti).idx     )
			end
			fx = SlicedSeparableSum(fxi,idxs)
		end
		fs = (fs...,fx)
	end
	if length(fs) > 1
		return SeparableSum(fs)  ##probably change constructor in Prox?
	else
		return fs[1]
	end
end

extract_proximable(xAll::Variable, t::Term) =  extract_functions(t)
extract_proximable(xAll::NTuple{N,Variable}, t::Term) where {N} = extract_proximable(xAll,(t,))

