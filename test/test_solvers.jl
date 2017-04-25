srand(0)

##############################################################################
### Regularized least squares, with two variable blocks to make things weird
println("Testing: regularized least squares, with two variable blocks to make things weird")

m, n1, n2 = 30, 50, 100

x1_fpg = Variable(n1)
x2_fpg = Variable(n2)

A1 = randn(m, n1)
A2 = randn(m, n2)
b = randn(m)

lam1 = 0.2
lam2 = 1.0

expr = ls(A1*x1_fpg + A2*x2_fpg - b) + lam1*norm(x1_fpg, 1) + lam2*norm(x2_fpg, 2)
prob = problem(expr)
@time sol = solve(prob, FPG(tol=1e-10,verbose=0))
println(sol)

x1_zerofpr = Variable(n1)
x2_zerofpr = Variable(n2)
expr = ls(A1*x1_zerofpr + A2*x2_zerofpr - b) + lam1*norm(x1_zerofpr, 1) + lam2*norm(x2_zerofpr, 2)
prob = problem(expr)
@time sol = solve(prob, ZeroFPR(tol=1e-10,verbose=0))
println(sol)

res = A1*~x1_fpg + A2*~x2_fpg - b
grad1 = A1'*res
grad2 = A2'*res
ind1_zero = (~x1_fpg .== 0)
subgr1 = lam1*sign(~x1_fpg)
subdiff1_low, subdiff1_upp = copy(subgr1), copy(subgr1)
subdiff1_low[ind1_zero] = -lam1
subdiff1_upp[ind1_zero] = +lam1
subgr2 = lam2*(~x2_fpg/norm(~x2_fpg, 2))

@test maximum(subdiff1_low + grad1) <= 1e-6
@test maximum(-subdiff1_upp - grad1) <= 1e-6
@test norm(grad2 + subgr2) <= 1e-6

##############################################################################
### Lasso problem with known solution
println("Testing: lasso problem with known solution")

m, n, nnz_x_star = 200, 100, 10
A = randn(m, n)
lam = 1.0
x_star = randn(n)
x_star[nnz_x_star+1:end] = 0.0
y_star = lam*sign(x_star)
b = A*x_star + A'\y_star
@test norm(A'*(A*x_star - b) + lam*sign(x_star)) <= 1e-12

# Solve with PG

x_pg = Variable(n)
expr = ls(A*x_pg - b) + lam*norm(x_pg, 1)
prob = problem(expr)
@time sol = solve(prob, PG(tol=1e-10,verbose=0))
println(sol)

@test norm(~x_pg - x_star, Inf) <= 1e-8
@test norm(A'*(A*~x_pg - b) + lam*sign(~x_pg)) <= 1e-6

# Solve with FPG

x_fpg = Variable(n)
expr = ls(A*x_fpg - b) + lam*norm(x_fpg, 1)
prob = problem(expr)
@time sol = solve(prob, FPG(tol=1e-10,verbose=0))
println(sol)

@test norm(~x_fpg - x_star, Inf) <= 1e-8
@test norm(A'*(A*~x_fpg - b) + lam*sign(~x_fpg)) <= 1e-6

# Solve with ZeroFPR

x_zerofpr = Variable(n)
expr = ls(A*x_zerofpr - b) + lam*norm(x_zerofpr, 1)
prob = problem(expr)
@time sol = solve(prob, ZeroFPR(tol=1e-10,verbose=0))
println(sol)

@test norm(~x_zerofpr - x_star, Inf) <= 1e-8
@test norm(A'*(A*~x_zerofpr - b) + lam*sign(~x_zerofpr)) <= 1e-5

################################################################################
### Problem with smooth, non-quadratic term
println("Testing: problem with smooth, non-quadratic term")

m, n, nnz_x_orig = 200, 500, 10
A = randn(m, n)
lam = 1.0
x_orig = randn(n)
x_orig[nnz_x_orig+1:end] = 0.0
b = A*x_orig + randn(m)

# Solve with PG

x_pg = Variable(n)
expr = smooth(norm(A*x_pg - b, 2)) + lam*norm(x_pg, 1)
prob = problem(expr)
@time sol = solve(prob, PG(tol=1e-8,verbose=0))
println(sol)

# Solve with FPG

x_fpg = Variable(n)
expr = smooth(norm(A*x_fpg - b, 2)) + lam*norm(x_fpg, 1)
prob = problem(expr)
@time sol = solve(prob, FPG(tol=1e-8,verbose=0))
println(sol)

# Solve with ZeroFPR

x_zerofpr = Variable(n)
expr = smooth(norm(A*x_zerofpr - b, 2)) + lam*norm(x_zerofpr, 1)
prob = problem(expr)
@time sol = solve(prob, ZeroFPR(tol=1e-8,verbose=0))
println(sol)
