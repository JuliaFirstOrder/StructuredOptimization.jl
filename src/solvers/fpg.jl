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


"""
# Fast Proximal Gradient Solver

## Usage


* `slv = FPG()` creates a `Solver` object that can be used in the function `solve`.
* Can be used with convex regularizers only.
* After solving a problem use `show(slv)` to visualize number of iterations, fixed point residual value, cost funtion value and time elapsed.


## Keyword Arguments

* `tol::Float64=1e-8`: tolerance
* `maxit::Int64=10000`: maximum number of iterations
* `verbose::Int64=1`: `0` verbose off, `1` print every 100 iteration, `2` print every iteration
* `stp_cr::Function=halt`: stopping criterion function
  * the function must have the following structure:

   `myhalt(slv::ForwardBackwardSolver,normfpr0::Float64,costprev::Float64)`

    * where `normfpr0` is the fixed point residual at x0
    * and `costprev` is the cost function value at the previous iteration
    * example: `myhalt(slv,normfpr0,FBE,FBEx) = slv.normfpr<slv.tol`

* `gamma::Float64=Inf`: stepsize γ, if γ = Inf upper bound is computed using:

  γ = || x0-(x0+ɛ) || / || ∇f(x0) - ∇f(x0+ɛ) ||

* `linesearch::Bool=true`: activates linesearch on stepsize γ

"""
FPG(; tol::Float64 = 1e-8,
      maxit::Int64 = 10000,
      verbose::Int64 = 1,
      stp_cr::Function = halt,
      linesearch::Bool = true,
      gamma::Float64 = Inf) =
FPG(tol, maxit, verbose, stp_cr, gamma,  0, Inf, Inf, NaN, linesearch, "Fast Proximal Gradient")

function solve!(L::Function, Ladj::Function, b::AbstractArray, g::ProximableFunction, x::AbstractArray, slv::FPG)

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
	xprev    = deepcopy(x)
	resxprev = deepcopy(resx)

	for slv.it = 1:slv.maxit

		# stopping criterion
		if slv.stp_cr(slv, normfpr0, costprev) break end

		# extrapolation
		y = x + slv.it/(slv.it+3) * (x - xprev)
		resy = resx + slv.it/(slv.it+3) * (resx - resxprev)

		# update iterates
		xprev =  deepcopy(x)
		resxprev = deepcopy(resx)
		costprev = copy(slv.cost)

		# compute gradient and f(y)
		fy = 0.5*vecnorm(resy)^2
		grady = Ladj(resy)

		# line search on gamma
		for j = 1:32
			gz = prox!(g, y - slv.gamma*grady, x, slv.gamma)
			fpr = y-x
			slv.normfpr = myVecnorm(fpr)
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

function solve(L::Function, Ladj::Function, b::AbstractArray, g::ProximableFunction, x0::AbstractArray, slv::FPG)
	x = deepcopy(x0) #copy initial conditions
	x, slv = solve!(L,Ladj,b,g,x,slv)
	return x, slv
end
