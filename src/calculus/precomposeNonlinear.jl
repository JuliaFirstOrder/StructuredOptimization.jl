import ProximalOperators: gradient! # this can be removed when moved to Prox

export PrecomposeNonlinear

struct PrecomposeNonlinear{P <: ProximableFunction, T <: AbstractOperator} <: ProximableFunction
	g::P
    G::T
end

PrecomposeNonlinear(g::P, G::T) where {P, T} = PrecomposeNonlinear{P, T}(g, G)

is_smooth(f::PrecomposeNonlinear) = is_smooth(f.g)

function (f::PrecomposeNonlinear)(x)
    return f.g(f.G*x)
end

function gradient!(y, f::PrecomposeNonlinear, x)
    # not the most efficient thing
    z, v = gradient(f.g, f.G*x)
    J = Jacobian(f.G, x)
    y = Ac_mul_B!(y, J, z)
    return v
end
