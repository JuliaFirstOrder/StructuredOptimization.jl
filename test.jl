#using PyCall
#pygui(:qt)
#using PyPlot   ## install PyPlot package by running 'Pkg.add("PyPlot")'
               ## on a Julia terminal

include("prox.jl")
include("ista.jl")
include("fista.jl")

srand(123)                # initialize random seed
m, n = 200, 1000          # problem dimensions
A = randn(m, n)           # create matrix
x_orig = sprandn(n, 1, 0.2)  # create sparse vector
sigma = 0.05              # noise std
y = A*x_orig + sigma*randn(m)  # output y

x_ls = A\y                # least squares solution
lambda_max = norm(A'y, Inf) # when lambda = lambda_max then solution is x = 0
lambda = 0.01*lambda_max   #
tol = 1e-6
maxit = 10000

@time x_ista, it =   ista(A, y, l1norm(lambda), zeros(n), maxit, tol, 1)
@time x_fista, it = fista(A, y, l1norm(lambda), zeros(n), maxit, tol, 1)


println("   ls error: $(norm(x_orig-x_ls   ))")
println(" ista error: $(norm(x_orig-x_ista ))")
println("fista error: $(norm(x_orig-x_fista))")



#figure()
#plot(x_orig, label = "orig")
#plot(x_ls, label = "ls")
#plot(x_ista, label = "ista")
#plot(x_fista, label = "fista")
#legend()
