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

function extract_proximable{N,M}(xAll::NTuple{N,Variable}, t::NTuple{M,Term})
	fs = ()
	for x in xAll
		tx = () #terms containing x
		for ti in t
			if x in variables(ti)
				tx = (tx...,ti)
			end
		end
		if isempty(tx)
			fx = IndFree()
		else
			fx = extract_proximable((x,),tx)
		end
		fs = (fs...,fx)
	end
	return SeparableSum(fs)
end

extract_proximable(xAll::Tuple{Variable}, t::Term) = extract_functions(t)
extract_proximable{N}(xAll::NTuple{N,Variable}, t::Term) = extract_proximable(xAll,(t,))
extract_proximable(xAll::Tuple{Variable}, t::Tuple{Term}) =  extract_proximable(xAll,t[1])

function extract_proximable{M}(xAll::Tuple{Variable}, t::NTuple{M,Term})
	#this case should happen only when all Terms in t are GetIndex
	fs, idxs = [], []
	for ti in t
		push!(fs,extract_functions(ti))
		push!(idxs,operator(ti).idx)
	end
	return SlicedSeparableSum(fs,idxs)
end
