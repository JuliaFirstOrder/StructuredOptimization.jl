immutable SmoothedPrimal
	s::CostFunction          #smooth cost function
	n::CostFunction          #nonsmooth cost function
	p::ProximableFunction
	gamma0::AbstractArray
end

SmoothedPrimal(s::CostFunction, n::CostFunction, p::ProximableFunction) =
SmoothedPrimal(s,n,p,logspace(1, -6, 10))

function solve(P::SmoothedPrimal, slv::Solver = default_slv())  
	slv0 = []
	for gamma in P.gamma0
		slv.verbose != 0 && println("solving for γ₀ = $(gamma) ")
		smoothed = smooth(P.n,gamma)
		slv0 = solve(P.s+smoothed, P.p, slv)
	end
	return slv0
end

function Base.show(io::IO, P::SmoothedPrimal)
	println("Smoothed Primal Problem")
	println()
	println("Smooth Cost Function:")
	println()
	show(P.s)
	println()
	println("Smoothed Cost Function:")
	println()
	show(P.n)
	println()
	println("Proximable operators:")
	println()
	show(P.p)
	println()
end
