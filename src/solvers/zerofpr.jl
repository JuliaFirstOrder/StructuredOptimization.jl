type ZeroFPR <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	lbfgs::LBFGS.Storage
	stp_cr::Function
	gamma::Float64
	it::Int
	normfpr::Float64
	cost::Float64
	time::Float64
	linesearch::Bool
	name::AbstractString
end

ZeroFPR(; tol::Float64 = 1e-8, maxit::Int64 = 10000,
          mem::Int64 = 10, verbose::Int64 = 1, 
	  stp_cr::Function = halt, linesearch::Bool = true, gamma::Float64 = Inf) =
ZeroFPR(tol, maxit, verbose, LBFGS.create(mem), stp_cr, gamma,
        0, Inf, Inf, NaN, linesearch, "ZeroFPR")

function solve!(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::ZeroFPR)

	tic();

	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2

	if slv.gamma == Inf #compute upper bound for Lipschitz constant using fd
		slv.gamma = get_gamma0(L,x,b,fx)
	end
	
	beta = 0.05
	sigma = beta/(4*slv.gamma)
	rbar_prev = zeros(x)

	H0, tau = 1., 1.
	d = zeros(x)

	# compute least squares residual and gradient
	gradx = Ladj(resx)
	xbar, gxbar = prox(g, x-slv.gamma*gradx, slv.gamma)
	r = x - xbar
	slv.normfpr = vecnorm(r)
	uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
	FBEx = uppbnd + gxbar

	fxbar, normfpr0, FBEprev, = NaN, NaN, NaN
	xbarbar   = copy(x)
	xbar_prev = copy(x)
	rbar_prev = copy(x)

	for slv.it = 1:slv.maxit

		# stopping criterion
		# TODO make this slv.stp_cr(slv::)
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, slv.normfpr, FBEprev, FBEx) break end
		FBEprev = copy(FBEx)

		resxbar = L(xbar) - b
		fxbar = 0.5*vecnorm(resxbar)^2

		# line search on gamma
		if slv.linesearch == true
			for j = 1:32
				if fxbar <= uppbnd break end
				slv.gamma = 0.5*slv.gamma
				sigma = 2*sigma
				gxbar = prox!(g, x-slv.gamma*gradx, slv.gamma, xbar)
				r = x - xbar
				slv.normfpr = vecnorm(r)
				resxbar = L(xbar) - b
				fxbar = 0.5*vecnorm(resxbar)^2
				uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
			end
		end

		if slv.it == 1 normfpr0 = copy(slv.normfpr) end

		slv.cost = fxbar+gxbar
		# print out stuff
		print_status(slv)

		# compute rbar
		gradxbar = Ladj(resxbar)
		prox!(g, xbar - slv.gamma*gradxbar, slv.gamma, xbarbar)
		rbar = xbar - xbarbar

		# compute direction according to L-BFGS
		if slv.it == 1
			d = -rbar
			LBFGS.reset(slv.lbfgs)
		else
			s = tau*d
			y = r - rbar_prev
			ys = real(vecdot(s,y))
			if ys > 0
				H0 = ys/real(vecdot(y,y))
				LBFGS.push(slv.lbfgs, s, y)
			end
			d = -LBFGS.matvec(slv.lbfgs, H0, rbar)
		end

		# store xbar and rbar for later use
		copy!(xbar_prev, xbar)
		copy!(rbar_prev, rbar)

		# line search on tau
		level = FBEx - sigma*slv.normfpr^2
		tau = 1.0
		Ad = L(d)
		ATAd = Ladj(Ad)
		for j = 1:32
			x = xbar_prev + tau*d
			resx = resxbar + tau*Ad
			fx = 0.5*vecnorm(resx)^2
			gradx = gradxbar + tau*ATAd
			gxbar = prox!(g, x - slv.gamma*gradx, slv.gamma, xbar)
			r = x - xbar
			slv.normfpr = vecnorm(r)
			uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
			if uppbnd + gxbar <= level break end
			tau = 0.5*tau
		end

		# evaluate FBE at x
		FBEx = uppbnd + gxbar
	end

	print_status(slv, 2*(slv.verbose>0))

	slv.time = toq();

	return xbar, slv
end

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x0::Array, slv::ZeroFPR)
	x = copy(x0) #copy initial conditions
	x, slv = solve!(L,Ladj,b,g,x,slv)
	return x, slv
end
