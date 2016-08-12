function fista(A::Array{Float64,2}, args...)

	L(y::Array{Float64}) = A*y
	Ladj(y::Array{Float64}) = A'*y
	fista(L, Ladj, args...)

end

function fista(A::SparseMatrixCSC, args...)

	L(y::Array{Float64}) = A*y
	Ladj(y::Array{Float64}) = A'*y
	fista(L, Ladj, args...)

end

function fista(A::Array{Complex{Float64},2}, args...)

	L(y::Array{Complex{Float64},1}) = A*y
	Ladj(y::Array{Complex{Float64},1}) = A'*y
	fista(L, Ladj, args...)

end

function fista(L::Function, Ladj::Function, b::Array, proxg::Function, x::Array, maxit=10000, tol=1e-5, verbose=1)

	gamma = 100.0
	z = xprev = x
	normfpr = Inf
	k = 0

	# compute least squares residual and f(x)
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2

	# initialize variables
	xprev = x
	resxprev = resx
	z = x
	resz = resx

	for k = 1:maxit

		# extrapolation
		y = x + k/(k+3) * (x - xprev)
		resy = resx + k/(k+3) * (resx - resxprev)

		# compute gradient and f(y)
		fy = 0.5*vecnorm(resy)^2
		grady = Ladj(resy)

		# line search on gamma
		for j = 1:32
			z, = proxg(y - gamma*grady, gamma)
			fpr = y-z
			normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fy - real(vecdot(grady,fpr)) + 1/(2*gamma)*normfpr^2
			if fz <= uppbnd; break; end
			gamma = 0.5*gamma
		end

		# stopping criterion
		if normfpr <= tol break end

		# print out stuff
		print_status(k, gamma, normfpr, verbose)

		# update iterates
		xprev = x
		resxprev = resx
		x = z
		resx = resz

	end

	print_status(k, gamma, normfpr, 2*(verbose>0))
	return z, k

end
