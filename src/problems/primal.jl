immutable Primal{A<:AffineOperator}
	Op::A
	g::ProximableFunction
end

function problem{T<:SmoothTerm}(fi::T, smooth::Array{OptTerm,1}, nonsmooth::Array{OptTerm,1})

	x = fi.A.x #extract variables
	g = get_prox(x, nonsmooth, smooth)
	 
	Primal(fi.A, g)
end

solve(P::Primal, args...) = solve(P.Op,P.g, args...)


