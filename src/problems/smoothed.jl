immutable SmoothedPrimal
	s::CompositeFunction          #smooth cost function
	n::CompositeFunction          #nonsmooth cost function
	p::ProximableFunction
	gamma0::AbstractArray
end

SmoothedPrimal(s::CompositeFunction, n::CompositeFunction, p::ProximableFunction) =
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
	println(io, "Smoothed Primal Problem")
	println(io)
	println(io, "Smooth Cost Function:")
	println(io)
	show(io, P.s)
	println(io)
	println(io, "Smoothed Cost Function:")
	println(io)
	show(io, P.n)
	println(io)
	println(io, "Proximable operators:")
	println(io)
	show(P.p)
end
