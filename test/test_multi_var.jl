println("Testing a problem with multiple variable blocks")
verb = 0
srand(123)
n,m = 10,100
X = zeros(n,m)
X[:,2], X[:,5] = randn(n),randn(n)
x = zeros(n*m)
x[10] = 20
y = randn(n*m)

A = OptVar(zeros(n*m))+reshape(OptVar(zeros(n,m)),n*m) 

L =x-> x[1]+x[2][:]
#Ladj =y-> [y,reshape(y,n,m)]
#
@test vecnorm(L([x,X])-A*[x,X]) <1e-12
prox_col = [NormL1(2.), NormL1(0.5)]

g = SeparableSum(prox_col)

x0 = [randn(size(x)),randn(size(X))]
b = A*[x,X]
x_pg, ~ = solve(A, b, g, x0, PG(verbose = verb))
@time x_pg, ~ = solve(A, b, g, x0, PG(verbose = verb))
x_fpg, ~ = solve(A, b, g, x0, FPG(verbose = verb))
@time x_fpg, ~ = solve(A, b, g, x0, FPG(verbose = verb))
x_zerofpr, ~ = solve(A, b, g, x0, ZeroFPR(verbose = verb))
@time x_zerofpr, ~ = solve(A, b, g, x0, ZeroFPR(verbose = verb))

@test vecnorm(x_pg-x_fpg)/(1+norm(x_fpg)) < 1e-5
@test vecnorm(x_zerofpr-x_pg)/(1+norm(x_zerofpr)) < 1e-5

n,m = 100,50
srand(123)
#total LS
A = randn(n,m)
x_star = full(sprandn(m,1,0.2))[:]
y_star = A*x_star
y_m = y_star+randn(n)

L =x-> A*x[1]-x[2]
Op = A*OptVar(zeros(m))-OptVar(zeros(n))

xx,yy = randn(size(x_star)),randn(size(y_star))
@test vecnorm(L([xx,yy])-Op*[xx,yy]) <1e-12

prox_col = [NormL1(10.), IndBox(y_m-2.,y_m+2.)]
g = SeparableSum(prox_col)

x0 = [zeros(x_star),zeros(y_star)]
@time x_z, ~ = solve(Op, zeros(n), g, x0, ZeroFPR(tol = 1e-9,verbose = verb))
