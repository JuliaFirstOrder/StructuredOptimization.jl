immutable Primal{A<:AffineOperator}
	Op::A
	p::ProximableFunction
end

function problem{T<:SmoothTerm}(g::T, smooth::Array{OptTerm,1}, nonsmooth::Array{OptTerm,1})

	x = g.A.x #extract variables
	p = get_prox(x, nonsmooth, smooth)
	 
	Primal(g.A, p)
end

solve(P::Primal, args...) = solve(P.Op,P.p, args...)


