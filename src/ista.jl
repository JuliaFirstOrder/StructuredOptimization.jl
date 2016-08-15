ista(A, args...) = ista(x -> A*x, y -> A'*y, args...)

function ista(L::Function, Ladj::Function, b::Array, proxg::Function, x::Array, maxit=10000, tol=1e-5, verbose=1)

	gamma = 100.0
	z = xprev = x
	normfpr = Inf
	k = 0

	# compute least squares residual and f(x)
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2

	# initialize variables
	z = x
	resz = resx
	fz = fx

	for k = 1:maxit

		# compute gradient
		gradx = Ladj(resx)

		# line search on gamma
		for j = 1:32
			z, = proxg(x - gamma*gradx, gamma)
			fpr = x-z
			normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fx - real(vecdot(gradx,fpr)) + 1/(2*gamma)*normfpr^2
			if fz <= uppbnd; break; end
			gamma = 0.5*gamma
		end

		# stopping criterion
		if normfpr <= tol break end

		# print out stuff
		print_status(k, gamma, normfpr, verbose)

		# update iterates
		x = z
		resx = resz
		fx = fz

	end

	print_status(k, gamma, normfpr, 2*(verbose>0))
	return z, k

end
