export ZeroFPR

type ZeroFPR <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	mem::Int64
	halt::Function
	gamma::Float64
	it::Int
	fpr::AbstractArray
	normfpr::Float64
	cost::Float64
	time::Float64
	linesearch::Bool
	cnt_matvec::Int
	cnt_grad_q::Int
	cnt_grad_s::Int
	cnt_prox::Int
end

fun_name(S::ZeroFPR) = "ZeroFPR"

"""
# ZeroFPR Solver

Default solver of RegLS.

## Usage

* `slv = ZeroFPR()` creates a `Solver` object that can be used in the function `solve`.

* Can be used with convex and noncovex regularizers.

* After solving a problem use `show(slv)` to visualize number of iterations, fixed point residual value, cost funtion value and time elapsed.


## Keyword Arguments

* `tol::Float64=1e-8`: tolerance used in stopping criterion
* `maxit::Int64=10000`: maximum number of iterations
* `mem::Int64=5`: L-BFGS memory
* `verbose::Int64=1`: `0` verbose off, `1` print every 100 iteration, `2` print every iteration
* `stp_cr::Function`: custom stopping criterion function
  * this function may be specified by the user and must have the following structure:

    `myhalt(slv::ForwardBackwardSolver,normfpr0::Float64,Fcurr::Float64,Fprev::Float64)`

    * `normfpr0` is the fixed point residual at x0
    * `Fcurr` is the objective value at the current iteration
    * `Fprev` is the objective value at the previous iteration
    * example: `myhalt(slv,normfpr0,FBE,FBEx) = slv.normfpr < tol`

* `gamma::Float64=Inf`: stepsize γ, if γ = Inf upper bound is computed using:

  γ0 = || x0-(x0+ɛ) || / || ∇f(x0) - ∇f(x0+ɛ) ||

* `linesearch::Bool=true`: activates linesearch on stepsize γ
"""

ZeroFPR(;
	tol::Float64 = 1e-8,
	maxit::Int64 = 10000,
	mem::Int64 = 5,
	verbose::Int64 = 1,
	halt::Function = halt_default,
	gamma::Float64 = Inf,
	linesearch::Bool = true) =
ZeroFPR(tol, maxit, verbose, mem, halt, gamma, 0, [], Inf, Inf, NaN, linesearch, 0, 0, 0, 0)

function solve(cf::CompositeFunction, solver::ZeroFPR)
	error("not yet implemented")
end

function apply(slv::ZeroFPR, x0::AbstractArray,
	f::ProximableFunction, L::LinearOperator, b::AbstractArray,
	g::ProximableFunction)
	if is_quadratic(s)
		return apply(slv, x0, f, L, b, NullFunction(), NullOperator(), [], g)
	else
		return apply(slv, x0, NullFunction(), NullOperator(), [], f, L, b, g)
	end
end

function apply(slv::ZeroFPR, x0::AbstractArray,
	s::ProximableFunction, L_s::LinearOperator, b_s::AbstractArray,
	q::ProximableFunction, L_q::LinearOperator, b_q::AbstractArray,
	g::ProximableFunction)

	tic()

	# q, s = split_Quadratic(f)
	x = deepcopy(x0)

	lbfgs = LBFGS(x, slv.mem)
	beta = 0.05

	if !isempty(q)
		res_q_x = A_mul_B(L_q, x) + b_q
		grad_q_res, q_x = gradient(q, res_q_x)
		grad_q_x = At_mul_B(L_q, grad_q_res)
	else
		res_q_x = []
		q_x = 0.0
		grad_q_x = 0.0
	end

	if !isempty(s)
		res_s_x = A_mul_B(L_s, x) + b_s
		grad_s_res, s_x = gradient(s, res_s_x)
		grad_s_x = At_mul_B(L_s, grad_s_res)
	else
		res_s_x = []
		s_x = 0.0
	  grad_s_x = 0.0
	end

	f_x = q_x + s_x
	grad_f_x = grad_q_x + grad_s_x

	if slv.gamma == Inf
		# compute upper bound for Lipschitz constant
		res_q_eps = A_mul_B(L_q, x+sqrt(eps())) + b_q
		res_s_eps = A_mul_B(L_s, x+sqrt(eps())) + b_s
		grad_q_res_eps, = gradient(q, res_q_eps)
		grad_s_res_eps, = gradient(s, res_s_eps)
		grad_f_eps = At_mul_B(L_q, grad_q_res_eps) + At_mul_B(L_s, grad_s_res_eps)
		Lf = deepvecnorm(grad_f_eps - grad_f_x)/(sqrt(eps()*deeplength(x)))
		slv.gamma = (1-beta)/Lf
	end

	sigma = beta/(4*slv.gamma)
	tau = 1.0

	gradstep = x - slv.gamma*grad_f_x
	xbar, g_xbar = prox(g, gradstep, slv.gamma)
	fpr_x = x - xbar
	slv.fpr = fpr_x
	slv.normfpr = deepvecnorm(fpr_x)
	uppbnd = f_x - real(deepvecdot(grad_f_x, fpr_x)) + 1/(2*slv.gamma)*slv.normfpr^2
	FBE_x = uppbnd + g_xbar

	# initialize variables

	normfpr0, FBE_prev, = NaN, NaN
	if !isempty(q)
		res_q_xbar = deepcopy(res_q_x)
		grad_q_xbar = deepcopy(grad_q_x)
		Aq_d = deepcopy(res_q_x)
		grad_q_Aq_d = deepcopy(grad_q_res)
		grad_q_d = deepcopy(grad_q_x)
	end
	if !isempty(s)
		res_s_xbar = deepcopy(res_s_x)
		grad_s_xbar = deepcopy(grad_s_x)
		As_d = deepcopy(res_s_x)
	end
	d = deepcopy(x)
	xbarbar = deepcopy(x)
	fpr_xbar = deepcopy(x)
	xbar_prev = deepcopy(x)
	fpr_xbar_prev = deepcopy(x)

	for slv.it = 1:slv.maxit

		# stopping criterion

		if slv.halt(slv) break end

		FBE_prev = FBE_x

		if !isempty(q)
			A_mul_B!(res_q_xbar, L_q, xbar)
			res_q_xbar .+= b_q
			q_xbar = q(res_q_xbar)
		else
			q_xbar = 0.0
		end

		if !isempty(s)
			A_mul_B!(res_s_xbar, L_s, xbar)
			res_s_xbar .+= b_s
			s_xbar = s(res_s_xbar)
		else
			s_xbar = 0.0
		end

		f_xbar = q_xbar + s_xbar

		# backtrack gamma

		if slv.linesearch == true
			for j = 1:32
				if f_xbar <= uppbnd + 1e-6*abs(f_xbar) break end
				slv.gamma = 0.5*slv.gamma
				sigma = 2*sigma
				gradstep .= x .- (*).(slv.gamma, grad_f_x)
				g_xbar = prox!(xbar, g, gradstep, slv.gamma)
				fpr_x .= (-).(x, xbar)
				slv.normfpr = deepvecnorm(fpr_x)
				if !isempty(q)
					A_mul_B!(res_q_xbar, L_q, xbar)
					res_q_xbar .+= b_q
					q_xbar = q(res_q_xbar)
				else
					q_xbar = 0.0
				end
				if !isempty(s)
					A_mul_B!(res_s_xbar, L_s, xbar)
					res_s_xbar .+= b_s
					s_xbar = s(res_s_xbar)
				else
					s_xbar = 0.0
				end
				f_xbar = q_xbar + s_xbar
				uppbnd = f_x - real(deepvecdot(grad_f_x, fpr_x)) + (1-beta)/(2*slv.gamma)*slv.normfpr^2
			end
		end

		if slv.it == 1 normfpr0 = copy(slv.normfpr) end

		# evaluate f+g and FBE at x

		FBE_x = uppbnd + g_xbar
		slv.cost = f_xbar + g_xbar

		# print out stuff

		print_status(slv)

		# compute rbar

		if !isempty(q)
			# println("res_q_xbar  : $(size(res_q_xbar))")
			# println("grad_q_res  : $(size(grad_q_res))")
			gradient!(grad_q_res, q, res_q_xbar)
			At_mul_B!(grad_q_xbar, L_q, grad_q_res)
		else
		  grad_q_xbar = 0.0
		end

		if !isempty(s)
			gradient!(grad_s_res, s, res_s_xbar)
			At_mul_B!(grad_s_xbar, L_s, grad_s_res)
		else
		  grad_s_xbar = 0.0
		end

		gradstep .= xbar
		gradstep .-= (*).(slv.gamma, grad_q_xbar)
		gradstep .-= (*).(slv.gamma, grad_s_xbar)

		prox!(xbarbar, g, gradstep, slv.gamma)

		fpr_xbar .=(-).(xbar, xbarbar)

		# compute direction according to L-BFGS

		if slv.it == 1
			deepcopy!(d, -fpr_xbar)
		else
			update!(lbfgs, xbar, xbar_prev, fpr_xbar, fpr_xbar_prev)
			A_mul_B!(d, lbfgs, fpr_xbar)
		end

		# store xbar and fpr_xbar for later use

		deepcopy!(fpr_xbar_prev, fpr_xbar)
		deepcopy!(xbar_prev, xbar)

		# backtrack tau

		level = FBE_x - sigma*slv.normfpr^2
		tau = 1.0

		if !isempty(q)
			A_mul_B!(Aq_d, L_q, d)
			gradient!(grad_q_Aq_d, q, Aq_d)
			At_mul_B!(grad_q_d, L_q, grad_q_Aq_d)
		end

		if !isempty(s)
			A_mul_B!(As_d, L_s, d)
		end

		for j = 1:32
			x .= (+).(xbar_prev, (*).(tau, d))
			if !isempty(q)
				res_q_x .= (+).(res_q_xbar, (*).(tau, Aq_d))
				q_x = q(res_q_x)
				grad_q_x .= (+).(grad_q_xbar, (*).(tau, grad_q_d))
			else
				q_x = 0.0
				grad_q_x = 0.0
			end
			if !isempty(s)
				res_s_x .= (+).(res_s_xbar, (*).(tau, As_d))
				s_x = gradient!(grad_s_res, s, res_s_x)
				At_mul_B!(grad_s_x, L_s, grad_s_res)
			else
				s_x = 0.0
				grad_s_x = 0.0
			end
			f_x = q_x + s_x
			grad_f_x .= (+).(grad_q_x, grad_s_x)
			gradstep .= x
			gradstep .-= (*).(slv.gamma, grad_f_x)
			g_xbar = prox!(xbar, g, gradstep, slv.gamma)
			fpr_x .= (-).(x, xbar)
			slv.normfpr = deepvecnorm(fpr_x)
			uppbnd = f_x - real(deepvecdot(grad_f_x, fpr_x)) + 1/(2*slv.gamma)*slv.normfpr^2
			if uppbnd + g_xbar <= level break end
			tau = 0.5*tau
		end

	end

	print_status(slv, 2*(slv.verbose>0))
	deepcopy!(x0, x)

	slv.time = toq()

	return slv

end
