export ZeroFPR

type ZeroFPR <: ForwardBackwardSolver
	tol::Float64
	maxit::Int64
	verbose::Int64
	mem::Int64
	halt::Function
	gamma::Float64
	it::Int
	normfpr::Float64
	cost::Float64
	time::Float64
	adaptive::Bool
	cnt_matvec::Int
	cnt_grad_q::Int
	cnt_grad_s::Int
	cnt_prox::Int
end

ZeroFPR(;
	tol::Float64 = 1e-8,
	maxit::Int64 = 10000,
	mem::Int64 = 5,
	verbose::Int64 = 1,
	halt::Function = halt_default,
	gamma::Float64 = Inf,
	adaptive::Bool = true) =
ZeroFPR(tol, maxit, verbose, mem, halt, gamma, 0, Inf, Inf, NaN, adaptive, 0, 0, 0, 0)

function apply!(slv::ZeroFPR, x0::T,
	f::ProximableFunction, L::AbstractOperator,
	g::ProximableFunction) where {T <: Union{AbstractArray, Tuple}}
	# if is_quadratic(f)
		# return apply!(slv, x0, Nullable{ProximableFunction}(), Nullable{AbstractOperator}(), Nullable(f), Nullable(L), g)
	# else
		return apply!(slv, x0, Nullable(f), Nullable(L), Nullable{ProximableFunction}(), Nullable{AbstractOperator}(), g)
	# end
end

function apply!(slv::ZeroFPR, x0::T,
	s::ProximableFunction, L_s::AbstractOperator,
	q::ProximableFunction, L_q::AbstractOperator,
	g::ProximableFunction) where {T <: Union{AbstractArray, Tuple}}
		return apply!(slv, x0, Nullable(s), Nullable(L_s), Nullable(q), Nullable(L_q), g)
end

################################################################################
################################################################################
# ZeroFPR algorithm
################################################################################
################################################################################

function apply!(slv::ZeroFPR, x0::T,
	Ns::Nullable{FS}, NL_s::Nullable{LS}, # smooth term
	Nq::Nullable{FQ}, NL_q::Nullable{LQ}, # quadratic term
	g::ProximableFunction) where {T <: Union{AbstractArray, Tuple}, FQ <: ProximableFunction, LQ <: AbstractOperator, FS <: ProximableFunction, LS <: AbstractOperator}

	tic()

	x = deepcopy(x0)
	grad_f_x = deepcopy(x0)
	gradstep = deepcopy(x0)
	fpr_x = deepcopy(x0)
	xbar = deepcopy(x0)

	lbfgs = LBFGS(x, slv.mem)
	beta = 0.05

	if !isnull(Ns)
		s = get(Ns)
		L_s = get(NL_s)
		grad_s_x = deepcopy(x0)
		res_s_x = L_s*x
		grad_s_res, s_x = gradient(s, res_s_x)
		Ac_mul_B!(grad_s_x, L_s, grad_s_res)
	else
		s_x = 0.0
		grad_s_x = 0.0
	end

	f_x = s_x
	grad_f_x = grad_s_x

	if slv.gamma == Inf
		# compute upper bound for Lipschitz constant
		if !isnull(Ns)
			res_s_eps = L_s*(x .+ sqrt(eps()))
			grad_s_res_eps, = gradient(s, res_s_eps)
			grad_s_eps = L_s'*grad_s_res_eps
		end
		grad_f_eps = grad_s_eps
		Lf = deepvecnorm(grad_f_eps .- grad_f_x)/(sqrt(eps()*deeplength(x)))
		slv.gamma = (1-beta)/Lf
	end

	sigma = beta/(4*slv.gamma)
	tau = 1.0

	# forward-backward steps from x

	deepaxpy!(gradstep, x, -slv.gamma, grad_f_x)
	g_xbar = prox!(xbar, g, gradstep, slv.gamma)
	deepaxpy!(fpr_x, x, -1.0, xbar)

	slv.normfpr = deepvecnorm(fpr_x)
	uppbnd = f_x - real(deepvecdot(grad_f_x, fpr_x)) + 1/(2*slv.gamma)*slv.normfpr^2
	FBE_x = uppbnd + g_xbar

	# initialize variables

	normfpr0, FBE_prev, = NaN, NaN
	if !isnull(Ns)
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

		if !isnull(Ns)
			A_mul_B!(res_s_xbar, L_s, xbar)
			s_xbar = s(res_s_xbar)
		else
			s_xbar = 0.0
		end

		f_xbar = s_xbar

		# backtrack gamma

		if slv.adaptive == true

			for j = 1:32

				if f_xbar <= uppbnd + 1e-6*abs(f_xbar) break end
				slv.gamma = slv.gamma/2
				sigma = 2*sigma

				# forward-backward steps from x

				deepaxpy!(gradstep, x, -slv.gamma, grad_f_x)
				g_xbar = prox!(xbar, g, gradstep, slv.gamma)
				deepaxpy!(fpr_x, x, -1.0, xbar)

				slv.normfpr = deepvecnorm(fpr_x)
				if !isnull(Ns)
					A_mul_B!(res_s_xbar, L_s, xbar)
					s_xbar = gradient!(grad_s_res, s, res_s_xbar)
				else
					s_xbar = 0.0
				end
				f_xbar = s_xbar
				uppbnd = f_x - real(deepvecdot(grad_f_x, fpr_x)) + (1-beta)/(2*slv.gamma)*slv.normfpr^2

			end

		end

		if slv.it == 1 normfpr0 = slv.normfpr end

		# evaluate FBE at x and f+g at xbar

		FBE_x = uppbnd + g_xbar
		slv.cost = f_xbar + g_xbar

		# print out stuff

		print_status(slv)

		# forward-backward steps from xbar

		if !isnull(Ns)
			Ac_mul_B!(grad_s_xbar, L_s, grad_s_res)
		else
			grad_s_xbar = 0.0
		end

		grad_f_xbar = grad_s_xbar

		deepaxpy!(gradstep, xbar, -slv.gamma, grad_s_xbar)
		prox!(xbarbar, g, gradstep, slv.gamma)

		# compute rbar

		deepaxpy!(fpr_xbar, xbar, -1.0, xbarbar)

		# compute direction according to L-BFGS

		if slv.it > 1
			update!(lbfgs, xbar, xbar_prev, fpr_xbar, fpr_xbar_prev)
		end
		A_mul_B!(d, lbfgs, fpr_xbar)

		# store xbar and fpr_xbar for later use

		deepcopy!(fpr_xbar_prev, fpr_xbar)
		deepcopy!(xbar_prev, xbar)

		# backtrack tau

		level = FBE_x - sigma*slv.normfpr^2
		tau = 1.0

		if !isnull(Ns)
			A_mul_B!(As_d, L_s, d)
		end

		for j = 1:32
			deepaxpy!(x, xbar_prev, tau, d)
			if !isnull(Ns)
				deepaxpy!(res_s_x, res_s_xbar, tau, As_d)
				s_x = gradient!(grad_s_res, s, res_s_x)
				Ac_mul_B!(grad_s_x, L_s, grad_s_res)
			else
				s_x = 0.0
				grad_s_x = 0.0
			end
			f_x = s_x
			grad_f_x = grad_s_x
			deepaxpy!(gradstep, x, -slv.gamma, grad_f_x)
			g_xbar = prox!(xbar, g, gradstep, slv.gamma)
			deepaxpy!(fpr_x, x, -1.0, xbar)
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
