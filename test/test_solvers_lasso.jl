srand(123)
verb = 0
solvers = [ZeroFPR(verbose = verb),FPG(verbose = verb),PG(verbose = verb)]

for slv in solvers

	####test lasso
	println("testing lasso (checking optimality condition) with solver "*RegLS.fun_name(slv))
	m, n = 10, 50
	A = randn(m, n)
	b = randn(m)
	x = Variable(n)
	lam_max = norm(A'*b,Inf)
	lam = 0.05*lam_max
	sol = minimize(ls(A*x-b) + lam*norm(x,1), slv)
	x = Variable(n)
	@time sol = minimize(ls(A*x-b) + lam*norm(x,1), slv)
	show(sol)
	# compute proximal-gradient step manually, check fixed-point residual
	gstep = ~x - sol.gamma*(A'*(A*(~x)-b))
	pgstep = sign(gstep).*max(0, abs.(gstep)-lam*sol.gamma)
	@test norm(pgstep - (~x), Inf) <= 1e-6
	# project subgr onto the subdifferential of the L1-norm at x
	subgr = -A'*(A*(~x)-b)
	subgr_proj = min(max(subgr, -lam), lam)
	subgr_proj[(~x) .< 0] = -lam
	subgr_proj[(~x) .> 0] = lam
	@test norm(subgr - subgr_proj, Inf) <= 1e-5

end
