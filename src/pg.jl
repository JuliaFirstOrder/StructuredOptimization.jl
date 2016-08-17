immutable PG <: Solver
	tol::Float64
	maxit::Int64
	verbose::Int64
end

PG(; tol::Float64 = 1e-8, maxit::Int64 = 10000, verbose::Int64 = 1) =
	PG(tol, maxit, verbose)

function solve(L::Function, Ladj::Function, b::Array, g::Function, x::Array, solver::PG)

	gamma = 100.0
	z = xprev = x
	normfpr = Inf
	k = 0

	# compute least squares residual and f(x)
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	fz = fx
	gz = Inf

	# initialize variables
	z = x
	resz = resx
	fz = fx

	for k = 1:solver.maxit

		# compute gradient
		gradx = Ladj(resx)

		# line search on gamma
		for j = 1:32
			z, gz = g(x - gamma*gradx, gamma)
			fpr = x-z
			normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fx - real(vecdot(gradx,fpr)) + 1/(2*gamma)*normfpr^2
			if fz <= uppbnd; break; end
			gamma = 0.5*gamma
		end

		# stopping criterion
		if normfpr <= solver.tol break end

		# print out stuff
		print_status(k, gamma, normfpr, fz+gz, solver.verbose)

		# update iterates
		x = z
		resx = resz
		fx = fz

	end

	print_status(k, gamma, normfpr, fz+gz, 2*(solver.verbose>0))
	return z, k

end
