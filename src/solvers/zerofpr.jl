immutable ZeroFPR <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	lbfgs::LBFGS.Storage
	stp_cr::Function
end

ZeroFPR(; tol::Float64 = 1e-8, maxit::Int64 = 10000,
          mem::Int64 = 10, verbose::Int64 = 1,stp_cr::Function = halt ) =
	ZeroFPR(tol, maxit, verbose, LBFGS.create(mem),stp_cr)

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, solver::ZeroFPR)

	gamma = 100.0
	beta = 0.05
	sigma = beta/(4*gamma)
	normr = Inf
	normr0 = Inf
	H0 = 1.0
	tau = 1.0
	d = zeros(x)
	rbar_prev = zeros(x)
	k = 0

	# compute least squares residual and gradient
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	gradx = Ladj(resx)
	xbar, gxbar = prox(g, gamma, x-gamma*gradx)
	fxbar = Inf
	FBEprev = Inf
	r = x - xbar
	normr = vecnorm(r)
	uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*gamma)*normr^2

	for k = 1:solver.maxit
		resxbar = L(xbar) - b
		fxbar = 0.5*vecnorm(resxbar)^2

		# line search on gamma
		for j = 1:32
			if fxbar <= uppbnd break end
			gamma = 0.5*gamma
			sigma = 2*sigma
			xbar, gxbar = prox(g, gamma, x-gamma*gradx)
			r = x - xbar
			normr = vecnorm(r)
			resxbar = L(xbar) - b
			fxbar = 0.5*vecnorm(resxbar)^2
			uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*gamma)*normr^2
		end

		if k == 1 normr0 = normr end

		# evaluate FBE at x
		FBEx = uppbnd + gxbar

		# stopping criterion
		if solver.stp_cr(solver.tol, gamma, normr0, normr, FBEprev, FBEx) break end
		FBEprev = FBEx

		# print out stuff
		print_status(k, gamma, normr, fxbar+gxbar, solver.verbose)

		# compute rbar
		gradxbar = Ladj(resxbar)
		xbarbar, = prox(g, gamma, xbar - gamma*gradxbar)
		rbar = xbar - xbarbar

		# compute direction according to L-BFGS
		if k == 1
			d = -rbar
			LBFGS.reset(solver.lbfgs)
		else
			s = tau*d
			y = r - rbar_prev
			ys = real(vecdot(s,y))
			if ys > 0
				H0 = ys/real(vecdot(y,y))
				LBFGS.push(solver.lbfgs, s, y)
			end
			d = -LBFGS.matvec(solver.lbfgs, H0, rbar)
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
			xbar, gxbar = prox(g, gamma, x - gamma*gradx)
			r = x - xbar
			normr = vecnorm(r)
			uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*gamma)*normr^2
			if uppbnd + gxbar <= level break end
			tau = 0.5*tau
		end
	end
	print_status(k, gamma, normr, fxbar+gxbar, 2*(solver.verbose>0))
	return xbar, k
end
