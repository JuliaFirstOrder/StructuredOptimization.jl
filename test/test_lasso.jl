using Base.Test

A = [1.0  2.0 -1.0 -1.0; -2.0 -1.0  0.0 -1.0; 3.0  0.0  4.0 -1.0; -4.0 -1.0 -3.0  1.0; 5.0  3.0  2.0  3.0]'
b = [1.0, 2.0, 3.0, 4.0]
m, n = size(A)
lambda = 0.1*norm(A'*b, Inf)
g = l1norm(lambda)
x0 = zeros(n)
x_star = [-3.877278911564627e-01, 0, 0, 2.174149659863943e-02, 6.168435374149660e-01]
tol = 1e-9
maxit = 10000
verbose = 1
tol_test = 1e-4

@time x_ista,  it =  ista(A, b, g, x0, maxit, tol, verbose)
@test norm(x_ista-x_star, Inf)/norm(x_star, Inf) <= tol_test

@time x_fista, it = fista(A, b, g, x0, maxit, tol, verbose)
@test norm(x_fista-x_star, Inf)/norm(x_star, Inf) <= tol_test

@time x_zerofpr, it = zerofpr(A, b, g, x0, maxit, tol, verbose)
@test norm(x_zerofpr-x_star, Inf)/norm(x_star, Inf) <= tol_test
