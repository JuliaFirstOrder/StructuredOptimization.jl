export PG, FPG

type PG <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	halt::Function
	gamma::Float64
	it::Int
	normfpr::Float64
	cost::Float64
	time::Float64
	adaptive::Bool
	fast::Bool
	cnt_matvec::Int
	cnt_prox::Int
end

PG(;
	tol::Float64 = 1e-8,
	maxit::Int64 = 10000,
	verbose::Int64 = 1,
	halt::Function = halt_default,
	adaptive::Bool = true,
	fast::Bool = false,
	gamma::Float64 = Inf) =
PG(tol, maxit, verbose, halt, gamma,  0, Inf, Inf, NaN, adaptive, fast, 0, 0)

# alias for fast = true
FPG(;
	tol::Float64 = 1e-8,
	maxit::Int64 = 10000,
	verbose::Int64 = 1,
	halt::Function = halt_default,
	adaptive::Bool = true,
	gamma::Float64 = Inf) =
PG(tol = tol, maxit = maxit, verbose = verbose, halt = halt, adaptive = adaptive, fast = true, gamma = gamma)

solver_name(slv::PG) = "$(slv.fast ? "Fast" : "") Proximal Gradient"

################################################################################
################################################################################
# solve! method
################################################################################
################################################################################

function solve!(terms::Tuple, solver::PG)
	xAll = extract_variables(terms)
	# Separate smooth and nonsmooth
	smooth, nonsmooth = split_smooth(terms)
	if is_proximable(nonsmooth)
		solver.verbose == true && println("--------------------Solving Primal---------------")
		f = extract_functions(smooth)
		L = extract_operators(xAll,smooth)
		g = extract_proximable(xAll,nonsmooth)
		apply!(solver, ~xAll, f, L, g)
		return solver
	end
	#strongly = [t for t in terms if is_strongly_convex(t) == true]
	#nonstrongly = [t for t in terms if is_strongly_convex(t) == false]
	#if false # TODO: here, a condition for "easily conjugable" should go
	#	# Solving the DUAL
        #		solver.verbose == true && println("-------------------- Solving Dual ---------------")
	#	return solver
	#end
	error("Sorry, I cannot solve this problem")
end

################################################################################
################################################################################
# Proximal gradient algorithm
################################################################################
################################################################################

function apply!(slv::PG, x0::T, f::ProximableFunction, L::AbstractOperator, g::ProximableFunction) where {T <: Union{AbstractArray, Tuple}}

	tic()

	x = deepcopy(x0)

	grad_f_x = deepcopy(x0)
	res_x = L*x
	grad_f_res, f_x = gradient(f, res_x)
	Ac_mul_B!(grad_f_x, L, grad_f_res)
	slv.cnt_matvec += 2
	g_x = Inf
	cost_xprev = Inf
	normfpr0 = Inf

	if slv.gamma == Inf
		# compute upper bound for Lipschitz constant
		grad_f_x_eps = deepcopy(x0)
		res_x_eps = L*(x .+ sqrt(eps()))
		grad_f_res_eps, = gradient(f, res_x_eps)
		Ac_mul_B!(grad_f_x_eps, L, grad_f_res_eps)
		slv.cnt_matvec += 2
		Lf = deepvecnorm(grad_f_x .- grad_f_x_eps)/(sqrt(eps()*deeplength(x)))
		slv.gamma = 1.0/Lf
	end

	# initialize variables

	fpr = deepcopy(x)
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
			deepaxpy!(gradstep, y, -slv.gamma, grad_f_y)
			g_x = prox!(x, g, gradstep, slv.gamma)
			slv.cnt_prox += 1
			deepaxpy!(fpr, y, -1.0, x)
			slv.normfpr = deepvecnorm(fpr)
			A_mul_B!(res_x, L, x)
			f_x = f(res_x)
			slv.cnt_matvec += 1
			if slv.adaptive == false break end
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
			# y = x + (it-1)/(it+2) * (x - xprev)
			deepaxpy!(y, x, (slv.it-1)/(slv.it+2), x)
			deepaxpy!(y, y, -(slv.it-1)/(slv.it+2), xprev)
			# res_y = res_x + (it-1)/(it+2) * (res_x - res_xprev)
			deepaxpy!(res_y, res_x, (slv.it-1)/(slv.it+2), res_x)
			deepaxpy!(res_y, res_y, -(slv.it-1)/(slv.it+2), res_xprev)
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
