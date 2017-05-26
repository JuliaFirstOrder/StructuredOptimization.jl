export PG, FPG

type PG <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	halt::Function
	gamma::Float64
	it::Int
	fpr::AbstractArray
	normfpr::Float64
	cost::Float64
	time::Float64
	linesearch::Bool
	fast::Bool
	cnt_matvec::Int
	cnt_prox::Int
end

fun_name(S::PG) = S.fast ? "Fast Proximal Gradient" : "Proximal Gradient"

"""
# Proximal Gradient Solver

## Usage


* `slv = PG()` creates a `Solver` object that can be used in the function `solve`.
* Can be used with convex regularizers only.
* After solving a problem use `show(slv)` to visualize number of iterations, fixed point residual value, cost funtion value and time elapsed.


## Keyword Arguments

* `tol::Float64=1e-8`: tolerance
* `maxit::Int64=10000`: maximum number of iterations
* `verbose::Int64=1`: `0` verbose off, `1` print every 100 iteration, `2` print every iteration
* `halt::Function`: custom stopping criterion function
  * this function may be specified by the user and must have the following structure:

    `myhalt(slv::ForwardBackwardSolver,normfpr0::Float64,Fcurr::Float64,Fprev::Float64)`

    * `normfpr0` is the fixed point residual at x0
    * `Fcurr` is the objective value at the current iteration
    * `Fprev` is the objective value at the previous iteration
    * example: `myhalt(slv,normfpr0,FBE,FBEx) = slv.normfpr < tol`

* `gamma::Float64=Inf`: stepsize γ, if γ = Inf upper bound is computed using:

  γ = || x0-(x0+ɛ) || / || ∇f(x0) - ∇f(x0+ɛ) ||

* `linesearch::Bool=true`: activates linesearch on stepsize γ
* `fast::Bool=true`: switches between proximal gradient and fast proximal gradient
"""

PG(;
	tol::Float64 = 1e-8,
	maxit::Int64 = 10000,
	verbose::Int64 = 1,
	halt::Function = halt_default,
	linesearch::Bool = true,
	fast::Bool = false,
	gamma::Float64 = Inf) =
PG(tol, maxit, verbose, halt, gamma,  0, [], Inf, Inf, NaN, linesearch, fast, 0, 0)

# alias for fast = true
FPG(;
	tol::Float64 = 1e-8,
	maxit::Int64 = 10000,
	verbose::Int64 = 1,
	halt::Function = halt_default,
	linesearch::Bool = true,
	gamma::Float64 = Inf) =
PG(tol = tol, maxit = maxit, verbose = verbose, halt = halt, linesearch = linesearch, fast = true, gamma = gamma)

function get_proximable_functions(terms::Vararg{Term})
	# loops here are probably better than one-liners
	fs = []
	for i in 1:length(terms)
		if length(t.A.Ls) != 1 && !is_gram_diagonal(t.A.Ls[1]) return [] end
		for j in (i+1):length(terms)
			if !isempty(intersect(variables(terms[i]), variables(terms[j])))
				return []
			end
		end
		if is_diagonal(t.A.Ls[1])
			D = get_gram_diagonal(t.A.Ls[1])
			f = PrecomposeDiagonal(t.f, D, t.A.b)
			append!(fs, f)
		else
			fwd = x -> t.A.Ls[1].L*x
			adj = y -> t.A.Ls[1].L'*y
			D = get_gram_diagonal(t.A.Ls[1])
			f = PrecomposeDiagonal(t.f, fwd, adj, D, t.A.b)
			append!(fs, f)
		end
	end
	return []
end

function solve(terms::Vector{Term}, solver::PG)
	# Separate smooth and nonsmooth
	smooth = [t for t in terms if is_smooth(t) == true]
	nonsmooth = [t for t in terms if is_smooth(t) == false]
	if is_proximable(nonsmooth...)
		println("Solving the PRIMAL")
		return solver
	end
	strongly = [t for t in terms if is_strongly_convex(t) == true]
	nonstrongly = [t for t in terms if is_strongly_convex(t) == false]
	if false # TODO: here, a condition for "easily conjugable" should go
		# Solving the DUAL
		println("Solving the DUAL")
		return solver
	end
	error("Sorry, I cannot solve this problem")
end

function apply{T <: Union{AbstractArray, Tuple}}(slv::PG, x0::T, f::ProximableFunction, L::LinearOperator, g::ProximableFunction)

	tic()

	x = deepcopy(x0)

	grad_f_x = deepcopy(x0)
	res_x = L*x # + b
	grad_f_res, f_x = gradient(f, res_x)
	Ac_mul_B!(grad_f_x, L, grad_f_res)
	slv.cnt_matvec += 2
	g_x = Inf
	cost_xprev = Inf
	normfpr0 = Inf

	if slv.gamma == Inf
		# compute upper bound for Lipschitz constant
		grad_f_x_eps = deepcopy(x0)
		res_x_eps = L*(x+sqrt(eps())) # + b
		grad_f_res_eps, = gradient(f, res_x_eps)
		Ac_mul_B!(grad_f_x_eps, L, grad_f_res_eps)
		slv.cnt_matvec += 2
		Lf = deepvecnorm(grad_f_x-grad_f_x_eps)/(sqrt(eps()*deeplength(x)))
		slv.gamma = 1.0/Lf
	end

	# initialize variables

	fpr = deepcopy(x)
	slv.fpr = fpr
	y = deepcopy(x)
	xprev = deepcopy(x)
	res_y = deepcopy(res_x)
	res_xprev = deepcopy(res_x)
	f_y = f_x
	grad_f_y = grad_f_x
	gradstep = deepcopy(x)

	for slv.it = 1:slv.maxit

		# line search on gamma
		for j = 1:32
			gradstep .= y .- slv.gamma.*grad_f_y
			g_x = prox!(x, g, gradstep, slv.gamma)
			slv.cnt_prox += 1
			fpr .= y .- x
			slv.normfpr = deepvecnorm(fpr)
			A_mul_B!(res_x, L, x)
			# res_x .+= b
			f_x = f(res_x)
			slv.cnt_matvec += 1
			if slv.linesearch == false break end
			uppbnd = f_y - real(deepvecdot(grad_f_y, fpr)) + (0.5/slv.gamma)*(slv.normfpr^2)
			if f_x <= uppbnd + 1e-6*abs(f_y) break end
			slv.gamma = 0.5*slv.gamma
		end

		slv.cost = f_x + g_x

		# print out stuff

		print_status(slv)

		# stopping criterion

		if slv.halt(slv) break end

		# extrapolation

		if slv.fast
			y .= x .+ (slv.it-1)/(slv.it+2) .* (x .- xprev)
			res_y .= res_x .+ (slv.it-1)/(slv.it+2) .* (res_x .- res_xprev)
		else
			# no need to copy, just move references around
			y = x
			res_y = res_x
		end

		# compute gradient and f(y)

		f_y = gradient!(grad_f_res, f, res_y)
		Ac_mul_B!(grad_f_y, L, grad_f_res)

		slv.cnt_matvec += 1

		# update iterates

		x, xprev = xprev, x
		res_x, res_xprev = res_xprev, res_x
		costprev = slv.cost

	end

	print_status(slv, 2*(slv.verbose>0))
	deepcopy!(x0, x)

	slv.time = toq()

	return slv

end
