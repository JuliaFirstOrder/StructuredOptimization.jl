A = [1.0  2.0 -1.0 -1.0; -2.0 -1.0  0.0 -1.0; 3.0  0.0  4.0 -1.0; -4.0 -1.0 -3.0  1.0; 5.0  3.0  2.0  3.0]'
b = [1.0, 2.0, 3.0, 4.0]
m, n = size(A)
lambda = norm(A'*b, Inf)
g  = NormL1(0.01*lambda)
g2 = NormL1(0.001*lambda)

x0 = zeros(n)
x_star = [-3.877278911564627e-01, 0, 0, 2.174149659863943e-02, 6.168435374149660e-01]
tol = 1e-8
maxit = 10000
verb = 0

@printf("Testing warm start \n")

solver = PG(verbose = verb, tol = tol)
@time x_pg,  ipg =  solve(A, b, g, x0, solver)
@time x_pg2, ipg2 =  solve(A, b, g2, x0, PG(verbose = verb, tol = tol))
@time x_pg2_ws, ipg2_ws =  solve(A, b, g2, x_pg, solver)
@test ipg2.it >= ipg2_ws.it

solver = FPG(verbose = verb, tol = tol)
@time x_fpg,  ifpg =  solve(A, b, g, x0, solver)
@time x_fpg2, ifpg2 =  solve(A, b, g2, x0, FPG(verbose = verb, tol = tol))
@time x_fpg2_ws, ifpg2_ws =  solve(A, b, g2, x_fpg, solver)
@test ifpg2.it >= ifpg2_ws.it

solver = ZeroFPR(verbose = verb, tol = tol)
@time x_z,  iz =  solve(A, b, g, x0, solver)
@time x_z2, iz2 =  solve(A, b, g2, x0, ZeroFPR(verbose = verb, tol = tol))
@time x_z2_ws, iz2_ws =  solve(A, b, g2, x_z, solver)
@test iz2.it >= iz2_ws.it
