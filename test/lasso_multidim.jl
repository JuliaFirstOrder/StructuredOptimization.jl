using RegLS
using ProximalOperators
using Base.Test

m, n = 5, 3 #sources, receivers
F = 2*pi*collect(1:20)  #frequencies
d = randn(n,m)  #delays

function MulDel(x::Array{Complex{Float64},2},d::Array{Float64,2},F::Array{Float64,1},n::Int64,m::Int64)
	y = zeros(Complex{Float64},n,length(F))
	for i = 1:n, j=1:m, f = 1:length(F)
		y[i,f] += x[j,f]*exp(-im*2*pi*d[i,j]*F[f])
	end
	return y
end


L    = x -> MulDel(x, d,F,n,m)
Ladj = x -> MulDel(x,-d',F,m,n)

x = randn(m,length(F))+im*randn(m,length(F))
y = randn(n,length(F))+im*randn(n,length(F))
norm(vecdot(L(x),y)-vecdot(x,Ladj(y))) #test adjoint operator works

x_star = full(sprandn(m,length(F),0.01)+im*sprandn(m,length(F),0.01))

b = L(x_star)
lambda = 0.006*vecnorm(Ladj(b), Inf)
g = NormL1(lambda)
x0 = zeros(Complex{Float64} ,m,length(F))
tol = 1e-8
maxit = 100000
verb = 0
tol_test = 1e-2

@printf("Solving a complex multidim lasso instance (m = %d, n = %d)\n", m, n)

x_ista,  slv =  solve(L, Ladj, b, g, x0, PG(verbose = verb, tol = tol, maxit = maxit))
@time x_ista,  slv =  solve(L, Ladj, b, g, x0, PG(verbose = verb, tol = tol, maxit = maxit))
@test slv.it < maxit
@test norm(x_ista-x_star, Inf)/norm(x_star, Inf) <= tol_test

x_fista, slv = solve(L, Ladj, b, g, x0, FPG(verbose = verb, tol = tol, maxit = maxit))
@time x_fista, slv = solve(L, Ladj, b, g, x0, FPG(verbose = verb, tol = tol, maxit = maxit))
@test slv.it < maxit
@test norm(x_fista-x_star, Inf)/norm(x_star, Inf) <= tol_test

x_zerofpr, slv = solve(L, Ladj, b, g, x0, ZeroFPR(verbose = verb, tol = tol, maxit = maxit))
@time x_zerofpr, slv = solve(L, Ladj, b, g, x0, ZeroFPR(verbose = verb, tol = tol, maxit = maxit))
@test slv.it < maxit
@test norm(x_zerofpr-x_star, Inf)/norm(x_star, Inf) <= tol_test
