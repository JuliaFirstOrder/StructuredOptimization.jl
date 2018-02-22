import ProximalOperators: gradient! # this can be removed when moved to Prox

export PrecomposeNonlinear

struct PrecomposeNonlinear{P <: ProximableFunction, 
                           T <: AbstractOperator,
                           D, C
                          } <: ProximableFunction
	g::P
    G::T
    bufD::D
    bufC::C
    bufC2::C
end

function PrecomposeNonlinear(g::P, G::T) where {P, T} 
    bufD  = blockzeros(domainType(G),  size(G,2)) 
    bufC  = blockzeros(codomainType(G),size(G,1)) 
    bufC2 = blockzeros(codomainType(G),size(G,1)) 
    PrecomposeNonlinear{P, T, typeof(bufD), typeof(bufC)}(g, G, bufD, bufC, bufC2)
end

is_smooth(f::PrecomposeNonlinear) = is_smooth(f.g)

function (f::PrecomposeNonlinear)(x)
    return f.g(f.G*x)
end

function gradient!(y::D, f::PrecomposeNonlinear{P,T,D,C}, x::D) where {P,T,D,C}
    A_mul_B!(f.bufC, f.G, x)
    v = gradient!(f.bufC2, f.g, f.bufC)
    J = Jacobian(f.G, x)
    y = Ac_mul_B!(y, J, f.bufC2)
    return v
end
