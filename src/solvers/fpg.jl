type FPG <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	stp_cr::Function
	resxprev::Array
	xprev::Array
	gamma::Float64
	it::Int
	normfpr::Float64
	cost::Float64
	time::Float64
	name::AbstractString
end

FPG(; tol::Float64 = 1e-8, maxit::Int64 = 10000, verbose::Int64 = 1,
      stp_cr::Function = halt, gamma::Float64 = Inf) =
FPG(tol, maxit, verbose, stp_cr, [], [], gamma,  0, Inf, Inf, NaN, "Fast Proximal Gradient")

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::FPG)

	tic();

	normfpr0 = Inf
	k = 0

	# compute least squares residual and f(x)
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	fz = fx
	gz = Inf
	costprev = fz + gz

	# initialize variables
	if slv.gamma == Inf
		slv.gamma = 100.0
	end
	if slv.xprev == []
		slv.xprev = x
		slv.resxprev = resx
	end
	z = x
	resz = resx

	for k = 1:slv.maxit
	slv.it = k

		# stopping criterion
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, slv.normfpr, costprev, slv.cost) break end
		costprev = copy(slv.cost)

		# extrapolation
		y = x + k/(k+3) * (x - slv.xprev)
		resy = resx + k/(k+3) * (resx - slv.resxprev)

		# compute gradient and f(y)
		fy = 0.5*vecnorm(resy)^2
		grady = Ladj(resy)

		# line search on gamma
		for j = 1:32
			z, gz = prox(g, y - slv.gamma*grady, slv.gamma)
			fpr = y-z
			slv.normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fy - real(vecdot(grady,fpr)) + 1/(2*slv.gamma)*slv.normfpr^2
			if fz <= uppbnd; break; end
			slv.gamma = 0.5*slv.gamma
		end

		if k == 1 normfpr0 = slv.normfpr end

		slv.cost = fz + gz

		# print out stuff
		print_status(slv)

		# update iterates
		slv.xprev = x
		slv.resxprev = resx
		x = z
		resx = resz

	end

	print_status(slv, 2*(slv.verbose>0))

	T = toq();

	return z, slv

end
