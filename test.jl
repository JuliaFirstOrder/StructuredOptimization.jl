# using PyCall
# pygui(:qt)
# using PyPlot   ## install PyPlot package by running 'Pkg.add("PyPlot")'
#                ## on a Julia terminal

include("prox.jl")
include("ista.jl")
include("fista.jl")

srand(123)                # initialize random seed
m, n = 200, 1500          # problem dimensions
A = randn(m, n)           # create matrix
x_orig = sprandn(n, 1, 0.2)  # create sparse vector
sigma = 0.01              # noise std
y = A*x_orig + sigma*randn(m)  # output y

x_ls = A\y                # least squares solution
lambda_max = norm(A'y, Inf) # when lambda = lambda_max then solution is x = 0
lambda = 0.05*lambda_max   #
tol = 1e-6
maxit = 10000

tic()
x_ista, it = ista(A, y, l1norm(lambda), zeros(n), maxit, tol, 1)
toc()

tic()
x_fista, it = fista(A, y, l1norm(lambda), zeros(n), maxit, tol, 1)
toc()

return # to suppress output of last expression?

# figure()
# plot(x, label = "original")
# plot(xh, label = "naive")
# plot(x1, label = "lasso")
# legend()
