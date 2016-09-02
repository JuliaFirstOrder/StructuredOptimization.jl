type PG <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	stp_cr::Function
	gamma::Float64
end

PG(; tol::Float64 = 1e-8, maxit::Int64 = 10000, verbose::Int64 = 1, stp_cr::Function = halt) =
PG(tol, maxit, verbose, stp_cr, Inf)

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::PG)

	tic();

	if slv.gamma == Inf
		slv.gamma = 100.0
	end
	normfpr = Inf
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

		# compute gradient
		gradx = Ladj(resx)

		# line search on gamma
		for j = 1:32
			z, gz = prox(g, x - slv.gamma*gradx, slv.gamma)
			fpr = x-z
			normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fx - real(vecdot(gradx,fpr)) + 1/(2*slv.gamma)*normfpr^2
			if fz <= uppbnd; break; end
			slv.gamma = 0.5*slv.gamma
		end

		if k == 1 normfpr0 = normfpr end

		cost = fz + gz

		# stopping criterion
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, normfpr, costprev, cost) break end
		costprev = cost

		# print out stuff
		print_status(k, slv.gamma, normfpr, cost, slv.verbose)

		# update iterates
		x = z
		resx = resz
		fx = fz

	end

	print_status(k, slv.gamma, normfpr, fz+gz, 2*(slv.verbose>0))

	T = toq();

	return z, BasicInfo(k, slv.gamma, normfpr, fz+gz, T);

end
