srand(123)
verb = 0
maxit = 50000
solvers = [ZeroFPR(verbose = verb, maxit = maxit),
	   FPG(    verbose = verb, maxit = maxit),
	   PG(     verbose = verb, maxit = maxit)]


m, n, l = 10, 50, 30
A1 = randn(m, n)
A2 = randn(m, l)
b = randn(m)

println("\n ------------------------------------------------------------- \n")
println("testing lasso single variable (checking optimality condition) \n")
for slv in solvers

	####test lasso single variable
	x = Variable(n)
	lam_max = norm(A1'*b,Inf)
	lam = 0.05*lam_max
	sol = minimize(ls(A1*x-b) + lam*norm(x,1), slv)
	x = Variable(n)
	@time sol = minimize(ls(A1*x-b) + lam*norm(x,1), slv)
	println(sol)
	# compute proximal-gradient step manually, check fixed-point residual
	gstep = ~x - sol.gamma*(A1'*(A1*(~x)-b))
	pgstep = sign(gstep).*max(0, abs.(gstep)-lam*sol.gamma)
	@test norm(pgstep - (~x), Inf) <= 1e-6
	# project subgr onto the subdifferential of the L1-norm at x
	subgr = -A1'*(A1*(~x)-b)
	subgr_proj = min(max(subgr, -lam), lam)
	subgr_proj[(~x) .< 0] = -lam
	subgr_proj[(~x) .> 0] = lam
	@test norm(subgr - subgr_proj, Inf) <= 1e-5

end

println("\n ------------------------------------------------------------- \n")
println("testing lasso double variable (checking optimality condition) \n")
for slv in solvers

	####test lasso double variable
	lam1,lam2 = 0.1, 0.05
	x,y = Variable(n), Variable(l)
	sol = minimize(ls(A1*x+A2*y-b) + lam1*norm(x,1)+ lam2*norm(y,1), slv)
	x,y = Variable(n), Variable(l)
	@time sol = minimize(ls(A1*x+A2*y-b) + lam1*norm(x,1)+ lam2*norm(y,1), slv)
	println(sol)
	# compute proximal-gradient step manually, check fixed-point residual
	gxstep = ~x - sol.gamma*(A1'*(A1*~x+A2*~y-b))
	pxgstep = sign(gxstep).*max(0, abs.(gxstep)-lam1*sol.gamma)
	@test norm(pxgstep - (~x), Inf) <= 1e-6
	gystep = ~y - sol.gamma*(A2'*(A1*~x+A2*~y-b))
	pygstep = sign(gystep).*max(0, abs.(gystep)-lam2*sol.gamma)
	@test norm(pygstep - (~y), Inf) <= 1e-6
	# project subgr onto the subdifferential of the L1-norm at x
	subgrx = -A1'*(A1*~x+A2*~y-b)
	subgr_projx = min(max(subgrx, -lam1), lam1)
	subgr_projx[(~x) .< 0] = -lam1
	subgr_projx[(~x) .> 0] = lam1
	# project subgr onto the subdifferential of the L1-norm at y
	@test norm(subgrx - subgr_projx, Inf) <= 1e-5
	subgry = -A2'*(A1*~x+A2*~y-b)
	subgr_projy = min(max(subgry, -lam2), lam2)
	subgr_projy[(~y) .< 0] = -lam2
	subgr_projy[(~y) .> 0] = lam2
	@test norm(subgry - subgr_projy, Inf) <= 1e-5

end

A1 = randn(m, n)+im*randn(m, n)
A2 = randn(m, l)+im*randn(m, l)
b = randn(m)+im*randn(m)

println("\n ------------------------------------------------------------- \n")
println("testing lasso single variable complex (checking optimality condition) \n")
for slv in solvers

	####test lasso single variable
	x = Variable(Complex{Float64},n)
	lam_max = norm(A1'*b,Inf)
	lam = 0.05*lam_max
	sol = minimize(ls(A1*x-b) + lam*norm(x,1), slv)
	x = Variable(Complex{Float64},n)
	@time sol = minimize(ls(A1*x-b) + lam*norm(x,1), slv)
	println(sol)
	# compute proximal-gradient step manually, check fixed-point residual
	gstep = ~x - sol.gamma*(A1'*(A1*(~x)-b))
	pgstep = sign(gstep).*max(0, abs.(gstep)-lam*sol.gamma)
	@test norm(pgstep - (~x), Inf) <= 1e-6
	# project subgr onto the subdifferential of the L1-norm at x
	subgr = -A1'*(A1*(~x)-b)
	subgr_proj = copy(subgr)
	subgr_proj[abs(~x) .> 0] = sign((~x)[abs(~x) .> 0])*lam
	@test norm(subgr - subgr_proj, Inf) <= 1e-5

end

println("\n ------------------------------------------------------------- \n")
println("testing lasso double variable complex (checking optimality condition) \n")
for slv in solvers

	####test lasso single variable
	lam1, lam2 = 0.5,0.7
	x = Variable(Complex{Float64},n)
	y = Variable(Complex{Float64},l)
	sol = minimize(ls(A1*x+A2*y-b) + lam1*norm(x,1) + lam2*norm(y,1), slv)
	x = Variable(Complex{Float64},n)
	y = Variable(Complex{Float64},l)
	@time sol = minimize(ls(A1*x+A2*y-b) + lam1*norm(x,1) + lam2*norm(y,1), slv)
	println(sol)
	# compute proximal-gradient step manually, check fixed-point residual
	gstep = ~x - sol.gamma*(A1'*(A1*~x+A2*~y-b))
	pgstep = sign(gstep).*max(0, abs.(gstep)-lam1*sol.gamma)
	@test norm(pgstep - (~x), Inf) <= 1e-6

	gstep = ~y - sol.gamma*(A2'*(A1*~x+A2*~y-b))
	pgstep = sign(gstep).*max(0, abs.(gstep)-lam2*sol.gamma)
	@test norm(pgstep - (~y), Inf) <= 1e-6

	# project subgr onto the subdifferential of the L1-norm at x
	subgr = -A1'*(A1*~x+A2*~y-b)
	subgr_proj = copy(subgr)
	subgr_proj[abs(~x) .> 0] = sign((~x)[abs(~x) .> 0])*lam1
	@test norm(subgr - subgr_proj, Inf) <= 1e-5

	# project subgr onto the subdifferential of the L1-norm at y
	subgr = -A2'*(A1*~x+A2*~y-b)
	subgr_proj = copy(subgr)
	subgr_proj[abs(~y) .> 0] = sign((~y)[abs(~y) .> 0])*lam2
	@test norm(subgr - subgr_proj, Inf) <= 1e-5

end
