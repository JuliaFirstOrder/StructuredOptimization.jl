Random.seed!(0)

################################################################################
### Regularized least squares, with two variable blocks to make things weird
################################################################################

println("Testing: regularized least squares, with two variable blocks to make things weird")

m, n1, n2 = 30, 50, 100

A1 = randn(m, n1)
A2 = randn(m, n2)
b = randn(m)

lam1 = 0.2
lam2 = 1.0

# Solve with FPG

x1_fpg = Variable(n1)
x2_fpg = Variable(n2)
expr = ls(A1*x1_fpg + A2*x2_fpg - b) + lam1*norm(x1_fpg, 1) + lam2*norm(x2_fpg, 2)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.FPG(tol=1e-10,verbose=0,maxit=20000))
println(sol)

# Solve with ZeroFPR

x1_zerofpr = Variable(n1)
x2_zerofpr = Variable(n2)
expr = ls(A1*x1_zerofpr + A2*x2_zerofpr - b) + lam1*norm(x1_zerofpr, 1) + lam2*norm(x2_zerofpr, 2)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.ZeroFPR(tol=1e-10,verbose=0))
println(sol)

# Solve with PANOC

x1_panoc = Variable(n1)
x2_panoc = Variable(n2)
expr = ls(A1*x1_panoc + A2*x2_panoc - b) + lam1*norm(x1_panoc, 1) + lam2*norm(x2_panoc, 2)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.PANOC(tol=1e-10,verbose=0))
println(sol)

# Solve with minimize, use default solver/options

x1 = Variable(n1)
x2 = Variable(n2)
@time sol = @minimize ls(A1*x1 + A2*x2 - b) + lam1*norm(x1, 1) + lam2*norm(x2, 2)
println(sol)

@test norm(~x1_fpg - ~x1_zerofpr, Inf)/(1+norm(~x1_zerofpr, Inf)) <= 1e-6
@test norm(~x2_fpg - ~x2_zerofpr, Inf)/(1+norm(~x2_zerofpr, Inf)) <= 1e-6
@test norm(~x1_fpg - ~x1_panoc, Inf)/(1+norm(~x1_panoc, Inf)) <= 1e-6
@test norm(~x2_fpg - ~x2_panoc, Inf)/(1+norm(~x2_panoc, Inf)) <= 1e-6
@test norm(~x1 - ~x1_zerofpr, Inf)/(1+norm(~x1_zerofpr, Inf)) <= 1e-3
@test norm(~x2 - ~x2_zerofpr, Inf)/(1+norm(~x2_zerofpr, Inf)) <= 1e-3

res = A1*~x1_fpg + A2*~x2_fpg - b
grad1 = A1'*res
grad2 = A2'*res
ind1_zero = (~x1_fpg .== 0)
subgr1 = lam1*sign.(~x1_fpg)
subdiff1_low, subdiff1_upp = copy(subgr1), copy(subgr1)
subdiff1_low[ind1_zero] .= -lam1
subdiff1_upp[ind1_zero] .= +lam1
subgr2 = lam2*(~x2_fpg/norm(~x2_fpg, 2))

@test maximum(subdiff1_low + grad1) <= 1e-6
@test maximum(-subdiff1_upp - grad1) <= 1e-6
@test norm(grad2 + subgr2) <= 1e-6

###############################################################################
## Lasso problem with known solution
###############################################################################

println("Testing: lasso problem with known solution")

m, n, nnz_x_star = 200, 100, 10
A = randn(m, n)
lam = 1.0
x_star = randn(n)
x_star[nnz_x_star+1:end] .= 0.0
y_star = lam*sign.(x_star)
b = A*x_star + A'\y_star
@test norm(A'*(A*x_star - b) + lam*sign.(x_star)) <= 1e-12

# Solve with PG

x_pg = Variable(n)
expr = ls(A*x_pg - b) + lam*norm(x_pg, 1)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.PG(tol=1e-10,verbose=0))
println(sol)

@test norm(~x_pg - x_star, Inf) <= 1e-8
@test norm(A'*(A*~x_pg - b) + lam*sign.(~x_pg)) <= 1e-6

# Solve with FPG

x_fpg = Variable(n)
expr = ls(A*x_fpg - b) + lam*norm(x_fpg, 1)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.FPG(tol=1e-10,verbose=0))
println(sol)

@test norm(~x_fpg - x_star, Inf) <= 1e-8
@test norm(A'*(A*~x_fpg - b) + lam*sign.(~x_fpg)) <= 1e-6

# Solve with ZeroFPR

x_zerofpr = Variable(n)
expr = ls(A*x_zerofpr - b) + lam*norm(x_zerofpr, 1)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.ZeroFPR(tol=1e-10,verbose=0))
println(sol)

@test norm(~x_zerofpr - x_star, Inf) <= 1e-8
@test norm(A'*(A*~x_zerofpr - b) + lam*sign.(~x_zerofpr)) <= 1e-5

# Solve with PANOC

x_panoc = Variable(n)
expr = ls(A*x_panoc - b) + lam*norm(x_panoc, 1)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.PANOC(tol=1e-10,verbose=0))
println(sol)

@test norm(~x_panoc - x_star, Inf) <= 1e-8
@test norm(A'*(A*~x_panoc - b) + lam*sign.(~x_panoc)) <= 1e-5

################################################################################
### Problem with smooth, non-quadratic term
################################################################################

println("Testing: problem with smooth, non-quadratic term")

m, n, nnz_x_orig = 200, 500, 10
A = randn(m, n)
lam = 1.0
x_orig = randn(n)
x_orig[nnz_x_orig+1:end] .= 0.0
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
@time sol = solve(prob, StructuredOptimization.ZeroFPR(tol=1e-8,verbose=0))
println(sol)

# Solve with PANOC

x_panoc = Variable(n)
expr = smooth(norm(A*x_panoc - b, 2)) + lam*norm(x_panoc, 1)
prob = problem(expr)
@time sol = solve(prob, StructuredOptimization.PANOC(tol=1e-8,verbose=0))
println(sol)

# Solve with minimize, default solver/options

x = Variable(n)
@time sol = @minimize smooth(norm(A*x - b, 2)) + lam*norm(x, 1)
println(sol)

@test norm(~x_fpg - ~x_zerofpr, Inf)/(1+norm(~x_zerofpr, Inf)) <= 1e-6
@test norm(~x_fpg - ~x_panoc, Inf)/(1+norm(~x_panoc, Inf)) <= 1e-6
@test norm(~x - ~x_zerofpr, Inf)/(1+norm(~x_zerofpr, Inf)) <= 1e-3

################################################################################
### Box-constrained least-squares
################################################################################

println("Testing: box-constrained least-squares")

m, n = 500, 200
A = randn(m, n)
lb, ub = -1.0, 1.0
x_orig = 2.0*randn(n)
x_orig = max.(lb, min.(ub, x_orig))
b = A*x_orig + randn(m)

# Solve with PG

x_pg = Variable(n)
expr = ls(A*x_pg - b)
prob = problem(expr, x_pg in [lb, ub])
@time sol = solve(prob, StructuredOptimization.PG(tol=1e-8,verbose=0))
println(sol)

@test norm(~x_pg - max.(lb, min.(ub, ~x_pg)), Inf) <= 1e-12
@test norm(~x_pg - max.(lb, min.(ub, ~x_pg - A'*(A*~x_pg - b))), Inf)/(1+norm(~x_pg, Inf)) <= 1e-8

# Solve with FPG

x_fpg = Variable(n)
expr = ls(A*x_fpg - b)
prob = problem(expr, x_fpg in [lb, ub])
@time sol = solve(prob, StructuredOptimization.FPG(tol=1e-8,verbose=0))
println(sol)

@test norm(~x_fpg - max.(lb, min.(ub, ~x_fpg)), Inf) <= 1e-12
@test norm(~x_fpg - max.(lb, min.(ub, ~x_fpg - A'*(A*~x_fpg - b))), Inf)/(1+norm(~x_fpg, Inf)) <= 1e-8

# Solve with ZeroFPR

x_zerofpr = Variable(n)
expr = ls(A*x_zerofpr - b)
prob = problem(expr, x_zerofpr in [lb, ub])
@time sol = solve(prob, StructuredOptimization.ZeroFPR(tol=1e-8,verbose=0))
println(sol)

@test norm(~x_zerofpr - max.(lb, min.(ub, ~x_zerofpr)), Inf) <= 1e-12
@test norm(~x_zerofpr - max.(lb, min.(ub, ~x_zerofpr - A'*(A*~x_zerofpr - b))), Inf)/(1+norm(~x_zerofpr, Inf)) <= 1e-8

# Solve with PANOC

x_panoc = Variable(n)
expr = ls(A*x_panoc - b)
prob = problem(expr, x_panoc in [lb, ub])
@time sol = solve(prob, StructuredOptimization.PANOC(tol=1e-8,verbose=0))
println(sol)

@test norm(~x_panoc - max.(lb, min.(ub, ~x_panoc)), Inf) <= 1e-12
@test norm(~x_panoc - max.(lb, min.(ub, ~x_panoc - A'*(A*~x_panoc - b))), Inf)/(1+norm(~x_panoc, Inf)) <= 1e-8

# Solve with minimize, default solver/options

x = Variable(n)
@time sol = @minimize ls(A*x - b) st x in [lb, ub]
println(sol)

@test norm(~x - max.(lb, min.(ub, ~x)), Inf) <= 1e-12
@test norm(~x - max.(lb, min.(ub, ~x - A'*(A*~x - b))), Inf)/(1+norm(~x, Inf)) <= 1e-4

################################################################################
### Non-negative least-squares from a known solution
################################################################################

println("Testing: non-negative least-squares from a known solution")

# Lagrangian:
#
#   0.5||Ax-b||^2 + y'(x-z) + [z >= 0]
#
# Optimality conditions:
#
#   A'(Ax-b) + y = 0, or A'b = A'Ax + y
#   x = z
#   z >= 0
#   y <= 0
#   y'z = 0

m, n, nnz_x_star = 500, 200, 100
A = randn(m, n)
x_star = rand(n)
x_star[nnz_x_star+1:end] .= 0.0
y_star = -rand(n)
y_star[1:nnz_x_star] .= 0.0
b = A*x_star + A'\y_star

# Solve with PG

x_pg = Variable(n)
expr = ls(A*x_pg - b)
prob = problem(expr, x_pg >= 0.0)
@time sol = solve(prob, StructuredOptimization.PG(tol=1e-8,verbose=0))
println(sol)

@test all(~x_pg .>= 0.0)
@test norm(~x_pg - x_star, Inf)/(1+norm(x_star, Inf)) <= 1e-8

# Solve with FPG

x_fpg = Variable(n)
expr = ls(A*x_fpg - b)
prob = problem(expr, x_fpg >= 0.0)
@time sol = solve(prob, StructuredOptimization.FPG(tol=1e-8,verbose=0))
println(sol)

@test all(~x_fpg .>= 0.0)
@test norm(~x_fpg - x_star, Inf)/(1+norm(x_star, Inf)) <= 1e-8

# Solve with ZeroFPR

x_zerofpr = Variable(n)
expr = ls(A*x_zerofpr - b)
prob = problem(expr, x_zerofpr >= 0.0)
@time sol = solve(prob, StructuredOptimization.ZeroFPR(tol=1e-8,verbose=0))
println(sol)

@test all(~x_zerofpr .>= 0.0)
@test norm(~x_zerofpr - x_star, Inf)/(1+norm(x_star, Inf)) <= 1e-8

# Solve with PANOC

x_panoc = Variable(n)
expr = ls(A*x_panoc - b)
prob = problem(expr, x_panoc >= 0.0)
@time sol = solve(prob, StructuredOptimization.PANOC(tol=1e-8,verbose=0))
println(sol)

@test all(~x_zerofpr .>= 0.0)
@test norm(~x_zerofpr - x_star, Inf)/(1+norm(x_star, Inf)) <= 1e-8

# Solve with minimize, default solver/options

x = Variable(n)
@time sol = @minimize ls(A*x - b) st x >= 0.0
println(sol)

@test all(~x .>= 0.0)
@test norm(~x - x_star, Inf)/(1+norm(x_star, Inf)) <= 1e-6
