type ZeroFPR <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	lbfgs::LBFGS.Storage
	stp_cr::Function
	sigma::Float64
	rbar_prev::Array
	gamma::Float64
	it::Int
	normfpr::Float64
	cost::Float64
	time::Float64
	name::AbstractString
end

ZeroFPR(; tol::Float64 = 1e-8, maxit::Int64 = 10000,
          mem::Int64 = 10, verbose::Int64 = 1, stp_cr::Function = halt, gamma::Float64 = Inf) =
ZeroFPR(tol, maxit, verbose, LBFGS.create(mem), stp_cr, 0., [], gamma, 
        0, Inf, Inf, NaN, "ZeroFPR")

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::ZeroFPR)

	tic();

	beta = 0.05

	if slv.gamma == Inf
		slv.gamma = 100.0
	end

	if slv.lbfgs.curridx == 0
		slv.sigma = beta/(4*slv.gamma)
		slv.rbar_prev = zeros(x)
	end
	H0, tau = 1., 1.
	d = zeros(x)

	# compute least squares residual and gradient
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	gradx = Ladj(resx)
	xbar, gxbar = prox(g, x-slv.gamma*gradx, slv.gamma)
	r = x - xbar
	slv.normfpr = vecnorm(r)
	uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
	FBEx = uppbnd + gxbar

	fxbar = normfpr0 = FBEprev = Inf

	for k = 1:slv.maxit
		slv.it = k

		# stopping criterion
		# TODO make this slv.stp_cr(slv::)
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, slv.normfpr, FBEprev, FBEx) break end

		resxbar = L(xbar) - b
		fxbar = 0.5*vecnorm(resxbar)^2

		# line search on gamma
		for j = 1:32
			if fxbar <= uppbnd break end
			slv.gamma = 0.5*slv.gamma
			slv.sigma = 2*slv.sigma
			xbar, gxbar = prox(g, x-slv.gamma*gradx, slv.gamma)
			r = x - xbar
			slv.normfpr = vecnorm(r)
			resxbar = L(xbar) - b
			fxbar = 0.5*vecnorm(resxbar)^2
			uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
		end

		if k == 1 normfpr0 = slv.normfpr end

		# evaluate FBE at x
		FBEx = uppbnd + gxbar

		FBEprev = FBEx

		# print out stuff
		slv.cost = fxbar+gxbar
		print_status(slv)

		# compute rbar
		gradxbar = Ladj(resxbar)
		xbarbar, = prox(g, xbar - slv.gamma*gradxbar, slv.gamma)
		rbar = xbar - xbarbar

		# compute direction according to L-BFGS
		if k == 1 && slv.lbfgs.curridx == 0
			d = -rbar
			LBFGS.reset(slv.lbfgs)
		else
			s = tau*d
			y = r - slv.rbar_prev
			ys = real(vecdot(s,y))
			if ys > 0
				H0 = ys/real(vecdot(y,y))
				LBFGS.push(slv.lbfgs, s, y)
			end
			d = -LBFGS.matvec(slv.lbfgs, H0, rbar)
		end

		# store xbar and rbar for later use
		xbar_prev = xbar
		slv.rbar_prev = rbar

		# line search on tau
		level = FBEx - slv.sigma*slv.normfpr^2
		tau = 1.0
		Ad = L(d)
		ATAd = Ladj(Ad)
		for j = 1:32
			x = xbar_prev + tau*d
			resx = resxbar + tau*Ad
			fx = 0.5*vecnorm(resx)^2
			gradx = gradxbar + tau*ATAd
			xbar, gxbar = prox(g, x - slv.gamma*gradx, slv.gamma)
			r = x - xbar
			slv.normfpr = vecnorm(r)
			uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
			if uppbnd + gxbar <= level break end
			tau = 0.5*tau
		end
	end

	print_status(slv, 2*(slv.verbose>0))

	slv.time = toq();

	return xbar, slv
end
