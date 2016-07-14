include("utils.jl")

function fista(A::Array{Float64,2}, args...)

	L(y::Array{Float64}) = A*y
	Ladj(y::Array{Float64}) = A'*y
	fista(L, Ladj, args...)

end

function fista(L::Function, Ladj::Function, b::Array{Float64}, proxg::Function, x::Array{Float64}, maxit=10000, tol=1e-5, verbose=1)

	gamma = 100.0
	z = xprev = x
	normfpr = Inf
	k = 0

	for k = 1:maxit

		# extrapolation
		y = x + k/(k+3) * (x - xprev)

		# compute least squares residual and gradient
		resy = L(y) - b
		fy = 0.5*vecnorm(resy)^2
		grady = Ladj(resy)

		# line search on gamma
		for j = 1:32
			gradstep = y - gamma*grady
			z, ~ = proxg(gradstep, gamma)
			fpr = y-z
			normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fy - dot(grady[:],fpr[:]) + 1/(2*gamma)*normfpr^2
			if fz <= uppbnd; break; end
			gamma = 0.5*gamma
		end

		# stopping criterion
		if normfpr <= tol break end

		# print out stuff
		print_status(k, gamma, normfpr, verbose)

		# update iterates
		xprev = x
		x = z

	end

	print_status(k, gamma, normfpr, 2*(verbose>0))
	return z, k

end
