using RegLS
using Base.Test

srand(123)
m, n = 10, 5 #sources, receivers
F = 2*pi*collect(1:30)  #frequencies
d = randn(n,m)  #delays

function MulDel(x::Array{Complex{Float64},2},d::Array{Float64,2},F::Array{Float64,1},n::Int64,m::Int64)
	y = zeros(Complex{Float64},length(F),n)
	for i = 1:n, j=1:m, f = 1:length(F)
		y[f,i] += x[f,j]*exp(-im*2*pi*d[i,j]*F[f])
	end
	return y
end

RealOrComplexArray = Union{Array{Float64}, Array{Complex{Float64}}}

L    = x -> MulDel(x, d,F,n,m)
Ladj = x -> MulDel(x,-d',F,m,n)

x = randn(length(F),m)+im*randn(length(F),m)
y = randn(length(F),n)+im*randn(length(F),n)
norm(vecdot(L(x),y)-vecdot(x,Ladj(y))) #test adjoint operator works

x_star = zeros(Complex{Float64},length(F),m)
x_star[:,randperm(m)[1:3]] = randn(length(F),3)+im*randn(length(F),3)

b = L(x_star)
lambda = 0.0001*vecnorm(Ladj(b), Inf)
g = indBallL20(3)
x0 = zeros(Complex{Float64} ,length(F),m)
tol = 1e-8
maxit = 100000
verbose = 1
tol_test = 1e-2

@printf("Solving a complex multidim lasso instance (m = %d, n = %d)\n", m, n)


@time x_ista,  it =  ista(L, Ladj, b, g, x0, maxit, tol, verbose)
@test it < maxit
@test norm(x_ista-x_star, Inf)/norm(x_star, Inf) <= tol_test

#@time x_fista, it = fista(L, Ladj, b, g, x0, maxit, tol, verbose)
#@test it < maxit
#@test norm(x_fista-x_star, Inf)/norm(x_star, Inf) <= tol_test

@time x_zerofpr, it = zerofpr(L, Ladj, b, g, x0, maxit, tol, verbose)
@test it < maxit
@test norm(x_zerofpr-x_star, Inf)/norm(x_star, Inf) <= tol_test
