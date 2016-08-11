using Base.Test

m, n = 300, 1000
A = sprandn(m, n, 0.1)
b = randn(m)
lambda = 0.05*norm(A'*b, Inf)
g = l1norm(lambda)
x0 = zeros(n)
tol = 1e-8
maxit = 100000
verbose = 1
tol_test = 1e-4

x_star, ~ = fista(A, b, g, x0, 100000, 1e-12, 0)

@time x_ista,  it =  ista(A, b, g, x0, maxit, tol, verbose)
@test norm(x_ista-x_star, Inf)/norm(x_star, Inf) <= tol_test

@time x_fista, it = fista(A, b, g, x0, maxit, tol, verbose)
@test norm(x_fista-x_star, Inf)/norm(x_star, Inf) <= tol_test

@time x_zerofpr, it = zerofpr(A, b, g, x0, maxit, tol, verbose)
@test norm(x_zerofpr-x_star, Inf)/norm(x_star, Inf) <= tol_test
