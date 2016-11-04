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

"""
# Proximal Gradient Solver 

## Usage 

* `slv = PG()` creates a `Solver` object that can be used in the function `solve`.
* Can be used with convex and noncovex regularizers.
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
	gradx = copy(Ladj(resx))
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

	for slv.it = 1:slv.maxit

		# stopping criterion
		if slv.stp_cr(slv, normfpr0, costprev) break end
		costprev = copy(slv.cost)

		# compute gradient
		gradx = Ladj(resx)

		# line search on gamma
		for j = 1:32
			gz = prox!(g, x - slv.gamma*gradx, z, slv.gamma)
			fpr = x-z
			slv.normfpr = vecnorm(fpr)
			resx = L(z) - b
			fz = 0.5*vecnorm(resx)^2
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
