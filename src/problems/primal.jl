immutable Primal
	s::CompositeFunction          #smooth cost function
	p::ProximableFunction
end

solve(P::Primal, slv::Solver = default_slv()) = solve(P.s, P.p, slv)

function Base.show(io::IO, P::Primal)
	println(io, "Primal Problem")
	println(io)
	println(io, "Smooth Cost Function:")
	println(io)
	show(io, P.s)
	println(io)
	println(io, "Proximable operators:")
	println(io)
	show(io, P.p)
end
