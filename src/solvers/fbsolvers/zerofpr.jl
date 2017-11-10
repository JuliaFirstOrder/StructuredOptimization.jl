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

solver_name(slv::ZeroFPR) = "ZeroFPR"

################################################################################
################################################################################
# solve! method
################################################################################
################################################################################

function solve!(terms::Tuple, solver::ZeroFPR)
	xAll = extract_variables(terms)
	# Separate smooth and nonsmooth
	smooth, nonsmooth = split_smooth(terms)
	# Separate quadratic and nonquadratic
	quadratic, smooth = split_quadratic(smooth)
	if is_proximable(nonsmooth)
		solver.verbose == true && println("--------------------Solving Primal---------------")
		g = extract_proximable(xAll,nonsmooth)
		if isempty(smooth)
			f = extract_functions(quadratic)
			L = extract_operators(xAll,quadratic)
			apply!(solver, ~xAll, f, L, g)
		elseif isempty(quadratic)
			f = extract_functions(smooth)
			L = extract_operators(xAll,smooth)
			apply!(solver, ~xAll, f, L, g)
		else
			f_q = extract_functions(quadratic)
			L_q = extract_operators(xAll,quadratic)
			f_s = extract_functions(smooth)
			L_s = extract_operators(xAll,smooth)
			apply!(solver, ~xAll, f_s, L_s, f_q, L_q, g)
		end
		return solver
	end
	error("Sorry, I cannot solve this problem")
end

################################################################################
################################################################################
# apply! shortcuts
################################################################################
################################################################################

function apply!(slv::ZeroFPR, x0::T,
	f::ProximableFunction, L::AbstractOperator,
	g::ProximableFunction) where {T <: Union{AbstractArray, Tuple}}
	if is_quadratic(f)
		return apply!(slv, x0, Nullable{ProximableFunction}(), Nullable{AbstractOperator}(), Nullable(f), Nullable(L), g)
	else
		return apply!(slv, x0, Nullable(f), Nullable(L), Nullable{ProximableFunction}(), Nullable{AbstractOperator}(), g)
	end
end

function apply!(slv::ZeroFPR, x0::T,
	s::ProximableFunction, Ls::AbstractOperator,
	q::ProximableFunction, Lq::AbstractOperator,
	g::ProximableFunction) where {T <: Union{AbstractArray, Tuple}}
		return apply!(slv, x0, Nullable(s), Nullable(Ls), Nullable(q), Nullable(Lq), g)
end

################################################################################
################################################################################
# ZeroFPR algorithm
################################################################################
################################################################################

function apply!(slv::ZeroFPR, x0::T,
	N_s::Nullable{FS}, N_Ls::Nullable{LS}, # smooth term
	N_q::Nullable{FQ}, N_Lq::Nullable{LQ}, # quadratic term
	g::ProximableFunction) where {T <: Union{AbstractArray, Tuple}, FQ <: ProximableFunction, LQ <: AbstractOperator, FS <: ProximableFunction, LS <: AbstractOperator}

	tic()

	x = deepcopy(x0)
	gradf_x = deepsimilar(x0)
	gradstep = deepsimilar(x0)
	fpr_x = deepsimilar(x0)
	xbar = deepsimilar(x0)

	lbfgs = LBFGS(x, slv.mem)
	beta = 0.05

	# compute gradient of f = s + q

	if !isnull(N_s)
		s = get(N_s)
		Ls = get(NL_s)
		grads_x = deepsimilar(x0)
		Ls_x = Ls*x
		grads_Ls_x, s_x = gradient(s, Ls_x)
		Ac_mul_B!(grads_x, Ls, grads_Ls_x)
	else
		s_x = 0.0
		grads_x = 0.0
	end
	if !isnull(N_q)
		q = get(N_q)
		Lq = get(N_Lq)
		gradq_x = deepsimilar(x0)
		Lq_x = Lq*x
		gradq_Lq_x, q_x = gradient(q, Lq_x)
		Ac_mul_B!(gradq_x, Lq, gradq_Lq_x)
		b, = gradient(q, deepzeros(Lq_x)) # TODO: REALLY NECESSARY?
	else
		q_x = 0.0
		gradq_x = 0.0
	end

	f_x = s_x + q_x
	gradf_x = grads_x .+ gradq_x

	# compute upper bound for Lipschitz constant (if needed)

	if slv.gamma == Inf
		xeps = x .+ sqrt(eps())
		if !isnull(N_s)
			Ls_xeps = Ls*xeps
			grads_Ls_xeps, = gradient(s, Ls_xeps)
			grads_xeps = Ls'*grads_Ls_xeps
		else
			grads_xeps = 0.0
		end
		if !isnull(N_q)
			Lq_xeps = Lq*xeps
			gradq_Lq_xeps, = gradient(q, Lq_xeps)
			gradq_xeps = Lq'*gradq_Lq_xeps
		else
			gradq_xeps = 0.0
		end
		gradf_xeps = grads_xeps .+ gradq_xeps
		Lf = deepvecnorm(gradf_xeps .- gradf_x)/(sqrt(eps()*deeplength(x)))
		slv.gamma = (1-beta)/Lf
	end

	sigma = beta/(4*slv.gamma)
	tau = 1.0

	# forward-backward step from x

	deepaxpy!(gradstep, x, -slv.gamma, gradf_x)
	g_xbar = prox!(xbar, g, gradstep, slv.gamma)
	deepaxpy!(fpr_x, x, -1.0, xbar)

	slv.normfpr = deepvecnorm(fpr_x)
	uppbnd = f_x - real(deepvecdot(gradf_x, fpr_x)) + 1/(2*slv.gamma)*slv.normfpr^2
	FBE_x = uppbnd + g_xbar

	# initialize variables

	normfpr0, FBE_prev = NaN, NaN
	if !isnull(N_s)
		Ls_xbar = deepsimilar(Ls_x)
		grads_Ls_xbar = deepsimilar(grads_Ls_x)
		grads_xbar = deepsimilar(grads_x)
		Ls_d = deepsimilar(Ls_x)
	end
	if !isnull(N_q)
		Lq_xbar = deepsimilar(Lq_x)
		gradq_Lq_xbar = deepsimilar(gradq_Lq_x)
		gradq_xbar = deepsimilar(gradq_x)
		Lq_d = deepsimilar(Lq_x)
		A_Lq_d = deepsimilar(gradq_Lq_x)
		Lqc_A_Lq_d = deepsimilar(gradq_x)
	end

	d = deepsimilar(x)
	xbarbar = deepsimilar(x)
	fpr_xbar = deepsimilar(x)
	xbar_prev = deepsimilar(x)
	fpr_xbar_prev = deepsimilar(x)
	f_xbar, q_xbar = Inf, Inf

	for slv.it = 1:slv.maxit

		# stopping criterion

		if slv.halt(slv) break end

		FBE_prev = FBE_x

		# backtrack gamma (if necessary)

		for j = 1:32

			if !isnull(N_s)
				A_mul_B!(Ls_xbar, Ls, xbar)
				s_xbar = gradient!(grads_Ls_xbar, s, Ls_xbar)
			else
				s_xbar = 0.0
			end
			if !isnull(N_q)
				A_mul_B!(Lq_xbar, Lq, xbar)
				q_xbar = gradient!(gradq_Lq_xbar, q, Lq_xbar)
			else
				q_xbar = 0.0
			end
			f_xbar = s_xbar + q_xbar
			uppbnd = f_x - real(deepvecdot(gradf_x, fpr_x)) + (1-beta)/(2*slv.gamma)*slv.normfpr^2

			if f_xbar <= uppbnd + 1e-6*abs(f_xbar) || slv.adaptive == false
				break
			end
			slv.gamma = slv.gamma/2
			sigma = 2*sigma

			# forward-backward steps from x

			deepaxpy!(gradstep, x, -slv.gamma, gradf_x)
			g_xbar = prox!(xbar, g, gradstep, slv.gamma)
			deepaxpy!(fpr_x, x, -1.0, xbar)

			slv.normfpr = deepvecnorm(fpr_x)

		end

		if slv.it == 1 normfpr0 = slv.normfpr end

		# evaluate FBE at x and f+g at xbar

		FBE_x = uppbnd + g_xbar
		slv.cost = f_xbar + g_xbar

		# print out stuff

		print_status(slv)

		# forward-backward steps from xbar

		if !isnull(N_s)
			Ac_mul_B!(grads_xbar, Ls, grads_Ls_xbar)
		else
			grads_xbar = 0.0
		end
		if !isnull(N_q)
			Ac_mul_B!(gradq_xbar, Lq, gradq_Lq_xbar)
		else
			gradq_xbar = 0.0
		end

		gradf_xbar = grads_xbar .+ gradq_xbar

		deepaxpy!(gradstep, xbar, -slv.gamma, gradf_xbar)
		prox!(xbarbar, g, gradstep, slv.gamma)

		# compute rbar

		deepaxpy!(fpr_xbar, xbar, -1.0, xbarbar)

		# compute direction according to L-BFGS

		if slv.it > 1
			update!(lbfgs, xbar, xbar_prev, fpr_xbar, fpr_xbar_prev)
		end
		A_mul_B!(d, lbfgs, fpr_xbar)

		# store xbar and fpr_xbar for later use

		fpr_xbar_prev, fpr_xbar = fpr_xbar, fpr_xbar_prev
		xbar_prev, xbar = xbar, xbar_prev

		# backtrack tau

		level = FBE_x - sigma*slv.normfpr^2
		tau = 1.0

		if !isnull(N_s)
			A_mul_B!(Ls_d, Ls, d)
		end
		if !isnull(N_q)
			A_mul_B!(Lq_d, Lq, d)
			gradient!(A_Lq_d, q, Lq_d)
			deepaxpy!(A_Lq_d, A_Lq_d, -1.0, b)
			Ac_mul_B!(Lqc_A_Lq_d, Lq, A_Lq_d)
			lin_coeff_q_d = real(deepvecdot(gradq_xbar, d))
			quad_coeff_q_d = real(deepvecdot(Lqc_A_Lq_d, d)/2)
		end

		for j = 1:32
			deepaxpy!(x, xbar_prev, tau, d)
			if !isnull(N_s)
				deepaxpy!(Ls_x, Ls_xbar, tau, Ls_d)
				s_x = gradient!(grads_Ls_x, s, Ls_x)
				Ac_mul_B!(grads_x, Ls, grads_Ls_x)
			else
				s_x = 0.0
				grads_x = 0.0
			end
			if !isnull(N_q)
				deepaxpy!(Lq_x, Lq_xbar, tau, Lq_d)
				deepaxpy!(gradq_Lq_x, gradq_Lq_xbar, tau, A_Lq_d)
				q_x = q_xbar + lin_coeff_q_d*tau + quad_coeff_q_d*(tau^2)
				deepaxpy!(gradq_x, gradq_xbar, tau, Lqc_A_Lq_d)
			else
				q_x = 0.0
				gradq_x = 0.0
			end
			f_x = s_x + q_x
			gradf_x = grads_x .+ gradq_x
			deepaxpy!(gradstep, x, -slv.gamma, gradf_x)
			g_xbar = prox!(xbar, g, gradstep, slv.gamma)
			deepaxpy!(fpr_x, x, -1.0, xbar)
			slv.normfpr = deepvecnorm(fpr_x)
			uppbnd = f_x - real(deepvecdot(gradf_x, fpr_x)) + 1/(2*slv.gamma)*slv.normfpr^2
			if uppbnd + g_xbar <= level break end
			tau = 0.5*tau
		end

	end

	print_status(slv, 2*(slv.verbose>0))
	deepcopy!(x0, x)

	slv.time = toq()

	return slv

end
