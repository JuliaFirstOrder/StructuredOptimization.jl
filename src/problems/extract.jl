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

# extract operators from terms

# returns all operators with an order dictated by xAll 

#single term, single variable
extract_operators(xAll::Tuple{Variable}, t::Term)  = operator(t)

extract_operators{N}(xAll::NTuple{N,Variable}, t::Term) = extract_operators(xAll, (t,))

#multiple terms, multiple variables
function extract_operators{N,M}(xAll::NTuple{N,Variable}, t::NTuple{M,Term})  
	ops = ()
	for ti in t
		xi   = variables(ti)
		opsi = operator(ti)
		ops = (ops..., sort_and_expand(xAll,xi,opsi))
	end
	return vcat(ops...)
end

function sort_and_expand{N}(xAll::NTuple{N,Variable}, xL::Tuple{Variable}, L::LinearOperator)
	ops = ()
	for i in eachindex(xAll)
		if xAll[i] == xL[1]
			ops = (ops...,L)
		else
			ops = (ops...,Zeros(eltype(~xAll[i]),size(xAll[i]),codomainType(L),size(L,1)))
		end
	end
	return hcat(ops...)
end

function sort_and_expand{N1,N2,M}(xAll::NTuple{N1,Variable}, xL::NTuple{N2,Variable}, L::HCAT{M,N2})
	ops = ()
	for i in eachindex(xAll)
		if xAll[i] in xL
			idx = findfirst(xAll[i].== xL)
			ops = (ops...,L[idx])
		else
			ops = (ops...,Zeros(eltype(~xAll[i]),size(xAll[i]),codomainType(L),size(L,1)))
		end
	end
	return HCAT(ops,L.mid,M)
end

# check terms are proximable
is_proximable{N}(xAll::NTuple{N,Variable}, t::Term) = is_gram_diagonal(t)

function is_proximable{N,M}(xAll::NTuple{N,Variable}, t::NTuple{M,Term})
	all(is_gram_diagonal(t)) == false && return false
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









