println("Testing a problem with multiple variable blocks")
verb = 0
#test with HCAT Op
srand(123)
n,m = 10,100
X = zeros(n,m)
X[:,2], X[:,5] = randn(n),randn(n)
x = zeros(n*m)
x[10] = 20
y = randn(n*m)

A = Variable(zeros(n*m))+reshape(Variable(zeros(n,m)),n*m) 
b = A*[x,X]
A = A-b

L =x-> x[1]+x[2][:]
#Ladj =y-> [y,reshape(y,n,m)]
#
@test vecnorm(L([x,X])-b-A*[x,X]) <1e-12
prox_col = [NormL1(2.), NormL1(0.5)]

g = SeparableSum(prox_col)

x0 = [randn(size(x)),randn(size(X))]
x_pg, ~ = solve(A, g, PG(verbose = verb))
@time x_pg, ~ = solve(A, g, PG(verbose = verb))
x_fpg, ~ = solve(A, g, FPG(verbose = verb))
@time x_fpg, ~ = solve(A, g, FPG(verbose = verb))
x_zerofpr, ~ = solve(A, g, ZeroFPR(verbose = verb))
@time x_zerofpr, ~ = solve(A, g, ZeroFPR(verbose = verb))

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
Op = A*Variable(zeros(m))-Variable(zeros(n))

xx,yy = randn(size(x_star)),randn(size(y_star))
@test vecnorm(L([xx,yy])-Op*[xx,yy]) <1e-12

prox_col = [NormL1(10.), IndBox(y_m-2.,y_m+2.)]
g = SeparableSum(prox_col)

x0 = [zeros(x_star),zeros(y_star)]
x_fpg, slv = solve(Op, g, FPG(tol = 1e-9,verbose = verb))
@time x_z, slv = solve(Op, g, FPG(tol = 1e-9,verbose = verb))
x_z, slv = solve(Op, g, ZeroFPR(tol = 1e-9,verbose = verb))
@time x_z, slv = solve(Op, g, ZeroFPR(tol = 1e-9,verbose = verb))

n,m = 100,50
srand(123)
# test with VCAT Op
X = Variable(randn(n,m))
b = [randn(n,m),randn(100)]

A = [dct(X);X[1:100]]-b
y = A*randn(n,m)

g = NormL1(0.5)

x_fpg, slv = solve(A, g, FPG(tol = 1e-9,verbose = verb))
@time x_fpg, slv = solve(A, g, FPG(tol = 1e-9,verbose = verb))
x_z, slv = solve(A, g, ZeroFPR(tol = 1e-9,verbose = verb))
@time x_z, slv = solve(A, g, ZeroFPR(tol = 1e-9,verbose = verb))
