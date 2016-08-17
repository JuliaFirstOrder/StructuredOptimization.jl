immutable FPG <: Solver
	tol::Float64
	maxit::Int64
	verbose::Int64
end

FPG(; tol::Float64 = 1e-8, maxit::Int64 = 10000, verbose::Int64 = 1) =
	FPG(tol, maxit, verbose)

function solve(L::Function, Ladj::Function, b::Array, g::Function, x::Array, solver::FPG)

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
	xprev = x
	resxprev = resx
	z = x
	resz = resx

	for k = 1:solver.maxit

		# extrapolation
		y = x + k/(k+3) * (x - xprev)
		resy = resx + k/(k+3) * (resx - resxprev)

		# compute gradient and f(y)
		fy = 0.5*vecnorm(resy)^2
		grady = Ladj(resy)

		# line search on gamma
		for j = 1:32
			z, gz = g(y - gamma*grady, gamma)
			fpr = y-z
			normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fy - real(vecdot(grady,fpr)) + 1/(2*gamma)*normfpr^2
			if fz <= uppbnd; break; end
			gamma = 0.5*gamma
		end

		# stopping criterion
		if normfpr <= solver.tol break end

		# print out stuff
		print_status(k, gamma, normfpr, fz+gz, solver.verbose)

		# update iterates
		xprev = x
		resxprev = resx
		x = z
		resx = resz

	end

	print_status(k, gamma, normfpr, fz+gz, 2*(solver.verbose>0))
	return z, k

end
