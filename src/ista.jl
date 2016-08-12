function ista(A::Array{Float64,2}, args...)

	L(y::Array{Float64}) = A*y
	Ladj(y::Array{Float64}) = A'*y
	ista(L, Ladj, args...)

end

function ista(A::SparseMatrixCSC, args...)

	L(y::Array{Float64}) = A*y
	Ladj(y::Array{Float64}) = A'*y
	ista(L, Ladj, args...)

end

function ista(L::Function, Ladj::Function, b::Array{Float64}, proxg::Function, x::Array{Float64}, maxit=10000, tol=1e-5, verbose=1)

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
			uppbnd = fx - vecdot(gradx,fpr) + 1/(2*gamma)*normfpr^2
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
