
using PyCall
pygui(:qt)
using PyPlot   ## install PyPlot package by running 'Pkg.add("PyPlot")' 
               ## on a Julia terminal
include("opt_fun.jl")

srand(123)           #initialize random seed
A = randn(50,100)    #create matrix
x = sprand(100,1,0.2)#create sparse vector x
y = A*x              #output y

xh = A\y             #naive estimation of x  
MaxIt,λ = 1000,5.    #num of iteration, λ parameter

x0 = zeros(length(x))#initialize x0 
TOL = 1e-5           #tolerance
x1,cf = fista(A,y,λ,x0,MaxIt,TOL,prox_l1) 

figure()
plot(x, label = "original")
plot(xh, label = "naive")
plot(x1, label = "lasso")
legend()

