# using PyCall
# pygui(:qt)
# using PyPlot   ## install PyPlot package by running 'Pkg.add("PyPlot")'
#                ## on a Julia terminal

include("fista.jl")

srand(123)                # initialize random seed
m, n = 50, 100            # problem dimensions
A = randn(m, n)           # create matrix
x_orig = sprandn(100, 1, 0.2)  # create sparse vector
sigma = 0.01              # noise std
y = A*x_orig + sigma*randn(m)  # output y

x_ls = A\y                # least squares solution
lambda_max = norm(A'y, Inf) # when lambda = lambda_max then solution is x = 0
lambda = 0.1*lambda_max   #

tic()
x_lasso, it = fista(A, y, lambda)
toc()

return # to suppress output of last expression?

# figure()
# plot(x, label = "original")
# plot(xh, label = "naive")
# plot(x1, label = "lasso")
# legend()
