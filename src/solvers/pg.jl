type PG <: ForwardBackwardSolver
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

PG(; tol::Float64 = 1e-8, 
     maxit::Int64 = 10000, 
     verbose::Int64 = 1,
     stp_cr::Function = halt, 
     gamma::Float64 = Inf, 
     linesearch::Bool = true) =
PG(tol, 
   maxit, 
   verbose, 
   stp_cr, 
   gamma, 
   0, Inf, Inf, NaN, linesearch, "Proximal Gradient")

function solve!(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x::Array, slv::PG)

	tic();

	# compute least squares residual and f(x)
	resx = L(x) - b
	gradx = Ladj(resx)
	fx = 0.5*vecnorm(resx)^2

	if slv.gamma == Inf #compute upper bound for Lipschitz constant using fd
		slv.gamma = get_gamma0(L,Ladj,x,gradx,b)
	end
	normfpr0 = Inf

	fz = copy(fx)
	gz = Inf
	costprev = Inf

	# initialize variables
	z = copy(x)
	resz = copy(resx)

	for slv.it = 1:slv.maxit

		# stopping criterion
		if slv.stp_cr(slv.tol, slv.gamma, normfpr0, slv.normfpr, costprev, slv.cost) break end
		costprev = copy(slv.cost)

		# compute gradient
		gradx = Ladj(resx)

		# line search on gamma
		for j = 1:32
			gz = prox!(g, x - slv.gamma*gradx, z, slv.gamma)
			fpr = x-z
			slv.normfpr = vecnorm(fpr)
			resz = L(z) - b
			fz = 0.5*vecnorm(resz)^2
			uppbnd = fx - real(vecdot(gradx,fpr)) + 1/(2*slv.gamma)*slv.normfpr^2
			if slv.linesearch == false; break; end
			if fz <= uppbnd; break; end
			slv.gamma = 0.5*slv.gamma
		end

		if slv.it == 1 normfpr0 = slv.normfpr end

		slv.cost = fz + gz

		# print out stuff
		print_status(slv)

		# update iterates
		copy!(x,z)
		copy!(resx,resz)
		fx = copy(fz)

	end

	print_status(slv, 2*(slv.verbose>0))

	slv.time = toq();

	return z, slv

end

function solve(L::Function, Ladj::Function, b::Array, g::ProximableFunction, x0::Array, slv::PG)
	x = copy(x0) #copy initial conditions
	x, slv = solve!(L,Ladj,b,g,x,slv)
	return x, slv
end
