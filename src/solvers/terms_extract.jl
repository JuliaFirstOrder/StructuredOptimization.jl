# returns all variables of a cost function, in terms of appearance
extract_variables(t::TermOrExpr) = variables(t) 

function extract_variables(t::NTuple{N,TermOrExpr}) where {N}
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
  f = displacement(t) == 0 ? t.f : PrecomposeDiagonal(t.f, 1.0, displacement(t)) #for now I keep this
  f = t.lambda == 1. ? f : Postcompose(f, t.lambda)                                  #for now I keep this
  #TODO change this
  return f
end
extract_functions(t::NTuple{N,Term}) where {N} = SeparableSum(extract_functions.(t))
extract_functions(t::Tuple{Term}) = extract_functions(t[1])

# extract functions from terms without displacement
function extract_functions_nodisp(t::Term)
  f = t.lambda == 1. ? t.f : Postcompose(t.f, t.lambda)
  return f
end
extract_functions_nodisp(t::NTuple{N,Term}) where {N} = SeparableSum(extract_functions_nodisp.(t))
extract_functions_nodisp(t::Tuple{Term}) = extract_functions_nodisp(t[1])

# extract operators from terms

# returns all operators with an order dictated by xAll

#single term, single variable
extract_operators(xAll::Tuple{Variable}, t::TermOrExpr)  = operator(t)
extract_operators(xAll::NTuple{N,Variable}, t::TermOrExpr) where {N} = extract_operators(xAll, (t,))

#multiple terms, multiple variables
function extract_operators(xAll::NTuple{N,Variable}, t::NTuple{M,TermOrExpr}) where {N,M}
  ops = ()
  for ti in t
    tex = expand(xAll,ti)
    ops = (ops...,sort_and_extract_operators(xAll,tex))
  end
  return vcat(ops...)
end

sort_and_extract_operators(xAll::Tuple{Variable}, t::TermOrExpr) = operator(t)

function sort_and_extract_operators(xAll::NTuple{N,Variable}, t::TermOrExpr) where {N}
  p = zeros(Int,N)
  xL = variables(t)
  for i in eachindex(xAll)
    p[i] = findfirst( xi -> xi == xAll[i], xL)
  end
  return operator(t)[p]
end

# extract affines from terms

# returns all affines with an order dictated by xAll

#single term, single variable
extract_affines(xAll::Tuple{Variable}, t::TermOrExpr)  = affine(t)

extract_affines(xAll::NTuple{N,Variable}, t::TermOrExpr) where {N} = extract_affines(xAll, (t,))

#multiple terms, multiple variables
function extract_affines(xAll::NTuple{N,Variable}, t::NTuple{M,TermOrExpr}) where {N,M}
  ops = ()
  for ti in t
    tex = expand(xAll,ti)
    ops = (ops...,sort_and_extract_affines(xAll,tex))
  end
  return vcat(ops...)
end

sort_and_extract_affines(xAll::Tuple{Variable}, t::TermOrExpr) = affine(t)

function sort_and_extract_affines(xAll::NTuple{N,Variable}, t::TermOrExpr) where {N}
  p = zeros(Int,N)
  xL = variables(t)
  for i in eachindex(xAll)
    p[i] = findfirst( xi -> xi == xAll[i], xL)
  end
  return affine(t)[p]
end

# expand term domain dimensions
function expand(xAll::NTuple{N,Variable}, t::Term) where {N}
  xt   = variables(t)
  C    = codomainType(operator(t))
  size_out = size(operator(t),1)
  ex = t.A

  for x in xAll
    if !( x in variables(t) ) 
      ex += Zeros(eltype(~x),size(x),C,size_out)*x
    end
  end
  return Term(t.lambda, t.f, ex)
end

function expand(xAll::NTuple{N,Variable}, ex::AbstractExpression) where {N}
  ex = convert(Expression,ex)
  xt   = variables(ex)
  C    = codomainType(operator(ex))
  size_out = size(operator(ex),1)

  for x in xAll
    if !( x in variables(ex) ) 
      ex += Zeros(eltype(~x),size(x),C,size_out)*x
    end
  end
  return ex
end

# extract function and merge operator
function extract_merge_functions(t::Term)
  if is_sliced(t)
    if typeof(operator(t)) <: Compose
      op = operator(t).A[2]
    else
      op = Eye(size(operator(t),1)...)
    end
  else
    op = operator(t)
  end
  if is_eye(op) 
    f = displacement(t) == 0 ? t.f : PrecomposeDiagonal(t.f, 1.0, displacement(t))
  elseif is_diagonal(op)
    f = PrecomposeDiagonal(t.f, diag(op), displacement(t))
  elseif is_AAc_diagonal(op)
    f = Precompose(t.f, op, diag_AAc(op), displacement(t))
  end
  f = t.lambda == 1. ? f : Postcompose(f, t.lambda) #for now I keep this
  #TODO change this
  return f
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
    else                            
      #multiple terms per variable
      #currently this happens only with GetIndex
      fxi,idxs = (),()
      for ti in tx
        fxi  = (fxi..., extract_merge_functions(ti))
        idx = typeof(operator(ti)) <: Compose ? operator(ti).A[1].idx : operator(ti).idx
        idxs = (idxs...,  idx   )
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

extract_proximable(xAll::Variable, t::Term) =  extract_merge_functions(t)
extract_proximable(xAll::NTuple{N,Variable}, t::Term) where {N} = extract_proximable(xAll,(t,))
