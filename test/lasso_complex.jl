using RegLS
using Prox
using Base.Test

srand(123)
m, n = 300, 500
A = randn(m, n)+im*randn(m,n)
b = randn(m)+im*randn(m)
lambda = 0.05*norm(A'*b, Inf)
g = NormL1(lambda)
x0 = zeros(Complex{Float64},n)
tol = 1e-8
maxit = 100000
verb = 0
tol_test = 1e-4

@printf("Solving a complex random lasso instance (m = %d, n = %d)\n", m, n)

x_star, ~ = solve(A, b, g, x0, FPG(verbose = verb))

@time x_ista,  slv =  solve(A, b, g, x0, PG(verbose = verb, tol = tol))
@test slv.it < maxit
@test norm(x_ista-x_star, Inf)/norm(x_star, Inf) <= tol_test

@time x_fista, slv = solve(A, b, g, x0, FPG(verbose = verb, tol = tol))
@test slv.it < maxit
@test norm(x_fista-x_star, Inf)/norm(x_star, Inf) <= tol_test

@time x_zerofpr, slv = solve(A, b, g, x0, ZeroFPR(verbose = verb, tol = tol))
@test slv.it < maxit
@test norm(x_zerofpr-x_star, Inf)/norm(x_star, Inf) <= tol_test
