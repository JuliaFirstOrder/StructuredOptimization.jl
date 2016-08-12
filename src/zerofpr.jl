function zerofpr(A::Array{Float64,2}, args...)

	L(y::Array{Float64}) = A*y
	Ladj(y::Array{Float64}) = A'*y
	zerofpr(L, Ladj, args...)

end

function zerofpr(A::SparseMatrixCSC, args...)

	L(y::Array{Float64}) = A*y
	Ladj(y::Array{Float64}) = A'*y
	zerofpr(L, Ladj, args...)

end

function zerofpr(L::Function, Ladj::Function, b::Array{Float64}, proxg::Function, x::Array{Float64}, maxit=10000, tol=1e-5, verbose=1)

	Lf = 1e-2
	beta = 0.05
	gamma = (1-beta)/Lf
	sigma = beta/(4*gamma)
	normr = Inf
	H0 = 1.0
	tau = 1.0
	d = zeros(size(x))
	rbar_prev = zeros(size(x))
	k = 0

	lbfgs = LBFGS.create(5)

	# compute least squares residual and gradient
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	gradx = Ladj(resx)
	xbar, gxbar = proxg(x-gamma*gradx, gamma)
	r = x - xbar
	normr = vecnorm(r)
	uppbnd = fx - vecdot(gradx,r) + 1/(2*gamma)*normr^2

	for k = 1:maxit
		resxbar = L(xbar) - b
		fxbar = 0.5*vecnorm(resxbar)^2

		# line search on gamma
		for j = 1:32
			if fxbar <= uppbnd break end
			Lf = 2*Lf
			gamma = 0.5*gamma
			sigma = 2*sigma
			xbar, gxbar = proxg(x-gamma*gradx, gamma)
			r = x - xbar
			normr = vecnorm(r)
			resxbar = L(xbar) - b
			fxbar = 0.5*vecnorm(resxbar)^2
			uppbnd = fx - vecdot(gradx,r) + 1/(2*gamma)*normr^2
		end

		# evaluate FBE at x
		FBEx = uppbnd + gxbar

		# stopping criterion
		if normr <= tol break end

		# print out stuff
		print_status(k, gamma, normr, verbose)

		# compute rbar
		gradxbar = Ladj(resxbar)
		xbarbar, gxbarbar = proxg(xbar - gamma*gradxbar, gamma)
		rbar = xbar - xbarbar

		# compute direction according to L-BFGS
		if k == 1
			d = -rbar
			LBFGS.reset(lbfgs)
		else
			s = tau*d
			y = r - rbar_prev
			ys = vecdot(s,y)
			if ys > 0
				H0 = ys/vecdot(y,y)
				LBFGS.push(lbfgs, s, y)
			end
			d = -LBFGS.matvec(lbfgs, H0, rbar)
		end

		# store xbar and rbar for later use
		xbar_prev = xbar
		rbar_prev = rbar

		# line search on tau
		level = FBEx - sigma*normr^2
		tau = 1.0
		Ad = L(d)
		ATAd = Ladj(Ad)
		for j = 1:32
			x = xbar_prev + tau*d
			resx = resxbar + tau*Ad
			fx = 0.5*vecnorm(resx)^2
			gradx = gradxbar + tau*ATAd
			xbar, gxbar = proxg(x - gamma*gradx, gamma)
			r = x - xbar
			normr = vecnorm(r)
			uppbnd = fx - vecdot(gradx,r) + 1/(2*gamma)*normr^2
			if uppbnd + gxbar <= level break end
			tau = 0.5*tau
		end
	end
	print_status(k, gamma, normr, 2*(verbose>0))
	return xbar, k
end
