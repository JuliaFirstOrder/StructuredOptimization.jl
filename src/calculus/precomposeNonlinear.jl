import ProximalOperators: gradient!, gradient # this can be removed when moved to Prox

export PrecomposeNonlinear

struct PrecomposeNonlinear{P <: ProximableFunction, 
                           T <: AbstractOperator,
                           D <: AbstractArray, 
                           C <: AbstractArray
                          } <: ProximableFunction
	g::P
  G::T
  bufD::D
  bufC::C
  bufC2::C
end

function PrecomposeNonlinear(g::P, G::T) where {P, T} 
  t, s = domainType(G), size(G,2)
  bufD = eltype(s) <: Int ? zeros(t,s) : ArrayPartition(zeros.(t,s))
  t, s = codomainType(G), size(G,1)
  bufC = eltype(s) <: Int ? zeros(t,s) : ArrayPartition(zeros.(t,s))
  bufC2 = eltype(s) <: Int ? zeros(t,s) : ArrayPartition(zeros.(t,s))
  PrecomposeNonlinear{P, T, typeof(bufD), typeof(bufC)}(g, G, bufD, bufC, bufC2)
end

is_smooth(f::PrecomposeNonlinear) = is_smooth(f.g)

function (f::PrecomposeNonlinear)(x)
    return f.g(f.G*x)
end

function gradient(f::PrecomposeNonlinear, x::ArrayPartition)
  y = zero(x)
  fy = gradient!(y,f,x)
  return y, fy
end

#TODO simplify this
function gradient!(y::D, f::PrecomposeNonlinear{P,T,D,C}, x::D) where {P,T,D <: ArrayPartition,C}
    mul!(f.bufC, f.G, x)
    v = gradient!(f.bufC2, f.g, f.bufC)
    J = Jacobian(f.G, x)
    y = mul!(y, J', f.bufC2)
    return v
end

function gradient!(y::D, f::PrecomposeNonlinear{P,T,D,C}, x::D) where {P,T,D <: AbstractArray,C}
    mul!(f.bufC, f.G, x)
    v = gradient!(f.bufC2, f.g, f.bufC)
    J = Jacobian(f.G, x)
    y = mul!(y, J', f.bufC2)
    return v
end
