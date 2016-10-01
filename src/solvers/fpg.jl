type FPG <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	stp_cr::Function
	gamma::Float64
	it::Int
	normfpr::Float64
	cost::Float64
	time::Float64
	linesearch::Bool
	name::AbstractString
end

FPG(; tol::Float64 = 1e-8, 
      maxit::Int64 = 10000, 
      verbose::Int64 = 1,
      stp_cr::Function = halt,
      linesearch::Bool = true,
      gamma::Float64 = Inf) =
FPG(tol, maxit, verbose, stp_cr, gamma,  0, Inf, Inf, NaN, linesearch, "Fast Proximal Gradient")

function solve!(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::FPG)

	tic();

	normfpr0 = NaN

	# compute least squares residual and f(x)
	resx = L(x) - b
	gradx = copy(Ladj(resx))
	fx = 0.5*vecnorm(resx)^2
	fz = copy(fx)
	gz = Inf
	costprev = fz + gz

	if slv.gamma == Inf #compute upper bound for Lipschitz constant using fd
		slv.gamma = get_gamma0(L,Ladj,x,gradx,b)
	end
		
	# initialize variables
	xprev    = copy(x)
	resxprev = copy(resx)

	for slv.it = 1:slv.maxit

		# stopping criterion
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, slv.normfpr, costprev, slv.cost) break end

		# extrapolation
		y = x + slv.it/(slv.it+3) * (x - xprev)
		resy = resx + slv.it/(slv.it+3) * (resx - resxprev)

		# update iterates
		copy!(xprev,      x)
		copy!(resxprev,resx)
		costprev = copy(slv.cost)

		# compute gradient and f(y)
		fy = 0.5*vecnorm(resy)^2
		grady = Ladj(resy)

		# line search on gamma
		for j = 1:32
			gz = prox!(g, y - slv.gamma*grady, x, slv.gamma)
			fpr = y-x
			slv.normfpr = vecnorm(fpr)
			resx = L(x) - b
			fz = 0.5*vecnorm(resx)^2
			uppbnd = fy - real(vecdot(grady,fpr)) + 1/(2*slv.gamma)*slv.normfpr^2
			if slv.linesearch == false; break; end
			if fz <= uppbnd; break; end
			slv.gamma = 0.5*slv.gamma
		end

		if slv.it == 1 normfpr0 = slv.normfpr end

		slv.cost = fz + gz

		# print out stuff
		print_status(slv)

	end

	print_status(slv, 2*(slv.verbose>0))

	slv.time = toq();

	return x, slv

end

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x0::Array, slv::FPG)
	x = copy(x0) #copy initial conditions
	x, slv = solve!(L,Ladj,b,g,x,slv)
	return x, slv
end
