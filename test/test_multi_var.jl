println("Testing a problem with multiple variable blocks")
verb = 0
srand(123)
n,m = 10,100
X = zeros(n,m)
X[:,2], X[:,5] = randn(n),randn(n)
x = zeros(n*m)
x[10] = 20
y = randn(n*m)

L =x-> x[1]+x[2][:]
Ladj =y-> [y,reshape(y,n,m)]

@test vecdot(L([x,X]),y)-vecdot([x,X],Ladj(y))<1e-12
prox_col = [NormL1(2.), NormL1(0.5)]

g = SeparableSum(prox_col)

x0 = [randn(size(x)),randn(size(X))]
b = L([x,X])
@time x_pg, ~ = solve(L,Ladj, b, g, x0, PG(verbose = verb))
@time x_fpg, ~ = solve(L,Ladj, b, g, x0, FPG(verbose = verb))
@time x_zerofpr, ~ = solve(L,Ladj, b, g, x0, ZeroFPR(verbose = verb))

@test vecnorm(x_pg-x_fpg)/(1+norm(x_fpg)) < 1e-5
@test vecnorm(x_zerofpr-x_pg)/(1+norm(x_zerofpr)) < 1e-5

n,m = 100,50
srand(123)
#total LS
A = randn(n,m)
x_star = full(sprandn(m,1,0.2))
y_star = A*x_star
y_m = y_star+randn(n)

L =x-> A*x[1]-x[2]
Ladj =y-> [A'*y,-y]

L2 =x-> A*x[1:m]-x[m+1:m+n]
Ladj2 =y-> [A'*y;-y]

xx,yy = randn(size(x_star)),randn(size(y_star))
zz = randn(n)

@test vecdot(L([xx,yy]),zz)-vecdot([xx,yy],Ladj(zz))<1e-12
@test vecdot(L2([xx;yy]),zz)-vecdot([xx;yy],Ladj2(zz))<1e-12
prox_col = [NormL1(10.), IndBox(y_m-2.,y_m+2.)]
g = SeparableSum(prox_col)
g2 = SlicedSeparableSum([prox_col[1]=>1:m,prox_col[2]=>m+1:m+n])

x0 = [zeros(x_star),zeros(y_star)]
@time x_z, ~ = solve(L, Ladj, zeros(n), g, x0, ZeroFPR(tol = 1e-9,verbose = verb))

x0 = [zeros(x_star);zeros(y_star)]
@time x_z2, ~ = solve(L2, Ladj2, zeros(n), g2, x0, ZeroFPR(tol = 1e-9,verbose = verb))

@test norm(x_z2-[x_z[1];x_z[2]])<1e-4


#x_fpg, ~ = solve(L, Ladj, zeros(n), g, x0, FPG(verbose = 1))

#using PyPlot
#figure()
#subplot(2,1,1)
#plot(x_z[1])
#plot(x_star)
#subplot(2,1,2)
#plot(x_z[2])
#plot(y_star)
