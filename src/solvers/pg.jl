type PG <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	stp_cr::Function
	gamma::Float64
	it::Int
	normfpr::Float64
	cost::Float64
	time::Float64
	name::AbstractString
end

PG(; tol::Float64 = 1e-8, maxit::Int64 = 10000, verbose::Int64 = 1, 
     stp_cr::Function = halt, gamma::Float64 = Inf) =
PG(tol, maxit, verbose, stp_cr, gamma, 0, Inf, Inf, NaN, "Proximal Gradient")

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::PG)

	tic();

	if slv.gamma == Inf
		slv.gamma = 100.0
	end
	normfpr0 = Inf
	k = 0

	# compute least squares residual and f(x)
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	fz = fx
	gz = Inf
	costprev = Inf

	# initialize variables
	z = x
	resz = resx

	for k = 1:slv.maxit
	slv.it = k

		# compute gradient
		gradx = Ladj(resx)

		# line search on gamma
		for j = 1:32
			z, gz = prox(g, x - slv.gamma*gradx, slv.gamma)
			fpr = x-z
			slv.normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fx - real(vecdot(gradx,fpr)) + 1/(2*slv.gamma)*slv.normfpr^2
			if fz <= uppbnd; break; end
			slv.gamma = 0.5*slv.gamma
		end

		if k == 1 normfpr0 = slv.normfpr end

		slv.cost = fz + gz

		# stopping criterion
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, slv.normfpr, costprev, slv.cost) break end
		costprev = copy(slv.cost)

		# print out stuff
		print_status(slv)

		# update iterates
		x = z
		resx = resz
		fx = fz

	end

	print_status(slv, 2*(slv.verbose>0))

	slv.time = toq();

	return z, slv

end
