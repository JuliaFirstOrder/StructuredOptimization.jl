include("utils.jl")

function fista(A, b, lambda)
	if typeof(A) == Array{Float64,2} fista(A, b, lambda, zeros(size(A, 2)))
	else error("A is not a matrix, provide the initial iterate x") end
end

function fista(A, b, lambda, x, maxit=10000, tol=1e-5, verbose=true)

	gamma = 100
	z = xprev = x
	normfpr = Inf
	k = 0

	for k = 1:maxit

		# extrapolation
		y = x + k/(k+3) * (x - xprev)

		# compute least squares residual and gradient
		resy = A*y - b
		grady = A'*resy

		# line search on gamma
		for j = 1:32
			gradstep = y - gamma*grady
			z = sign(gradstep).*max(0, abs(gradstep)-gamma*lambda)
			fpr = y-z
			normfpr = norm(fpr)
			resz = A*z - b
			fz = 0.5*norm(resz)^2
			uppbnd = 0.5*norm(resy)^2 - (grady'fpr)[1,1] + 1/(2*gamma)*normfpr^2
			if fz <= uppbnd; break; end
			gamma = 0.5*gamma
		end

		# print out stuff
		if verbose
			if k == 1 || k%1000 == 0 print_status() end
			if k == 1 || k%100 == 0 print_status(k, gamma, normfpr) end
		end

		# stopping criterion
		if normfpr <= tol break end

		# update iterates
		xprev = x
		x = z

	end

	if verbose print_status(k, gamma, normfpr) end
	return z, k

end
