m, n = 20, 50
r = 5
s = 0.1

B_orig = randn(m, r)*randn(r, n) + 0.1*randn(m, n)
M = rand(m, n);
M = 1.0*(M .<= s)
B = M.*B_orig
L = Ladj = X -> M.*X
g = IndBallRank(r)
x0 = randn(m, n)
A = diagop(Variable(x0),M)-B
verb = 0
maxit = 10000
tol = 1e-6

#println("$(norm(L(x0)-A*x0))")

@printf("Solving matrix completion (m = %d, n = %d)\n", m, n)

@time x_ista, slv = solve(A, g, PG(verbose = verb, tol = tol))
@test slv.it < slv.maxit
println(slv)

@time x_zerofpr, slv = solve(A, g, ZeroFPR(verbose = verb, tol = tol))
@test slv.it < slv.maxit
println(slv)
