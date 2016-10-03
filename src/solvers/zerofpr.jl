type ZeroFPR <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	mem::Int64
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
          mem::Int64 = 5, verbose::Int64 = 1, 
	  stp_cr::Function = halt, linesearch::Bool = true, gamma::Float64 = Inf) =
ZeroFPR(tol, maxit, verbose, mem, stp_cr, gamma,
        0, Inf, Inf, NaN, linesearch, "ZeroFPR")

function solve!(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::ZeroFPR)

	tic();
	lbfgs = LBFGS.create(slv.mem, x)

	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	gradx = copy(Ladj(resx))

	if slv.gamma == Inf #compute upper bound for Lipschitz constant using fd
		slv.gamma = get_gamma0(L,Ladj,x,gradx,b)
	end
	
	beta = 0.05
	sigma = beta/(4*slv.gamma)

	tau = 1.

	# compute least squares residual and gradient
	xbar, gxbar = prox(g, x-slv.gamma*gradx, slv.gamma)
	r = x - xbar
	slv.normfpr = vecnorm(r)
	uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
	FBEx = uppbnd + gxbar

	# initialize variables
	fxbar, normfpr0, FBEprev, = NaN, NaN, NaN
	rbar_prev = zeros(x)
	xbar_prev = zeros(x)
	xbarbar = x

	for slv.it = 1:slv.maxit

		# stopping criterion
		# TODO make this slv.stp_cr(slv::)
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, slv.normfpr, FBEprev, FBEx) break end
		FBEprev = copy(FBEx)

		resxbar = L(xbar) - b
		fxbar = 0.5*vecnorm(resxbar)^2

		# line search on gamma
		for j = 1:32
			if fxbar <= uppbnd break end
			slv.gamma = 0.5*slv.gamma
			sigma = 2*sigma
			gxbar = prox!(g, x-slv.gamma*gradx, xbar, slv.gamma)
			r[:] = x - xbar
			slv.normfpr = vecnorm(r)
			resxbar[:] = L(xbar) - b
			fxbar = 0.5*vecnorm(resxbar)^2
			uppbnd = fx - real(vecdot(gradx,r)) + 1/(2*slv.gamma)*slv.normfpr^2
		end

		if slv.it == 1 normfpr0 = copy(slv.normfpr) end

		# evaluate FBE at x
		FBEx = uppbnd + gxbar

		slv.cost = fxbar+gxbar
		# print out stuff
		print_status(slv)

		# compute rbar
		gradxbar = copy(Ladj(resxbar))
		prox!(g, xbar - slv.gamma*gradxbar, xbarbar, slv.gamma)
		rbar = xbar - xbarbar

		# compute direction according to L-BFGS
		if slv.it == 1
			LBFGS.push!(lbfgs, rbar)
		else
			LBFGS.push!(lbfgs, xbar, xbar_prev, rbar, rbar_prev)
		end

		# store xbar and rbar for later use
		copy!(rbar_prev, rbar)
		copy!(xbar_prev, xbar)

		# line search on tau
		level = FBEx - sigma*slv.normfpr^2
		tau = 1.0
		Ad = L(lbfgs.d)
		ATAd = Ladj(Ad)
		for j = 1:32
			x[:] = xbar_prev + tau*lbfgs.d
			resx[:] = resxbar + tau*Ad
			fx = 0.5*vecnorm(resx)^2
			gradx[:] = gradxbar + tau*ATAd
			gxbar = prox!(g, x - slv.gamma*gradx, xbar, slv.gamma)
			r[:] = x - xbar
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

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x0::Array, slv::ZeroFPR)
	x = copy(x0) #copy initial conditions
	x, slv = solve!(L,Ladj,b,g,x,slv)
	return x, slv
end
