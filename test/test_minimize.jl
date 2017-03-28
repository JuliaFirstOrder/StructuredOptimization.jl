#verb = 1
#slv = ZeroFPR(verbose = verb)
println("testing cost function constructors")

x,y,z,u = OptVar(randn(4)), OptVar(randn(4)), OptVar(randn(2,2)), OptVar(randn(4))
cf = ls(x)+ls(y)
@test length(variable(cf)) == 2
cf = ls(x+y)+ls(y)
@test length(variable(cf)) == 2
cf = ls(x+y)+ls(y)+ls(z)
@test length(variable(cf)) == 3
cf = ls(x)+ls(y)+ls(z)
@test length(variable(cf)) == 3
cf = 10*(4*ls(x)+ls(x+y)+2*ls(z))
@test length(variable(cf)) == 3
@test  cf.f[1].lambda == 40 && cf.f[2].lambda == 10 && cf.f[3].lambda == 20
show(cf)
cf = ls(x+y-randn(4))+ls(y+reshape(z,4)+u)+1*norm(randn(4).*u,1)+2*norm(z-randn(2,2),2)+3*norm(y,1)
@test length(variable(cf)) == 4
show(cf)
	
x_sorted = sort(variable(cf),by = object_id)

println("\n test split \n")
smooth, proximable, nonsmooth = RegLS.split(cf)
println("\n Smooth:")
show(smooth)
println("\n Proximable:")
show(proximable)
println("\n NonSmooth:")
show(nonsmooth)


println("\n test merge prox \n")
p = RegLS.mergeProx(x_sorted, proximable)
show(x_sorted)
show(p.fs)

println("\n test merge prox \n")
smooth_exp = RegLS.mergeSmooth(x_sorted, smooth)
show(RegLS.affine(smooth))
show(RegLS.affine(smooth_exp))

out1 = RegLS.affine(smooth_exp)[1](optData(x_sorted))
out2 = RegLS.affine(smooth)[1]([x.x,y.x])
@test norm(out1-out2) <= 1e-8

out1 = RegLS.affine(smooth_exp)[2](optData(x_sorted))
out2 = RegLS.affine(smooth)[2]([y.x,z.x,u.x])
@test norm(out1-out2) <= 1e-8
#norm(RegLS.affine(smooth_exp)[1]([X,Y,Z,U]))


#	
#
#
######test single block of variable
#println("testing single variable primal")
#x = OptVar(zeros(5))
#minimize(ls(fft(x)-fft(randn(5)))+1e-2*norm(x,1), slv)
#
#n,m = 5,3
#A = randn(m,n)
#x = OptVar(zeros(n))
#b1,b2 = randn(m),randn(n)
#minimize(ls(A*x-randn(m))+1e-2*norm(x,1), slv)
#
#x = OptVar(zeros(n))
#minimize(ls(x-b2)+1e-2*norm(x,1), slv)
#
#X, = minimize(ls(A*x-b1), [norm(5.0*x,1) <= 1/1e-2], slv)
#@test norm(X,1) <= 1/1e-2
#
#b = [b1,(A\b1)[1:2]]
#X, = minimize(ls([A*x; x[1:2]]-b), [norm(x,1) <= 1/1e-8], slv)
#@test norm(X,1) <= 1/1e-8
#@test vecnorm([A*X,X[1:2]]-b) <= 1e-5
#
#####test merge regularize
#X, = minimize(ls(A*x-b1)+1e7*ls(x-b2), [norm(x,1) <= 1/1e-2], slv)
#@test norm(X-b2) < 1e-5
#
#
#println("testing block variables primal")
#x, y = OptVar(n), OptVar(Complex{Float64},n)
#X, = minimize(ls(fft(x)+y-fft(b2))+1e-7*norm(x,1), slv)
#@test norm(fft(X[1])+X[2]-fft(b2))<1e-7
#
#x, y = OptVar(n), OptVar(m)
#X, = minimize(ls(A*x+y)+1e-2*norm(x,1), norm(y-b1,1) <= 1e-2, slv)
#
#@test (norm(X[2]-b1,1)-1e-2)<1e-8
#
#X, = minimize(ls(A*x+y), [norm(x,2)<=1e2, norm(2*y-b1,1) <= 1e-2], slv)
#
#@test (norm(2.*X[2]-b1,1)-1e-2)<1e-8
#@test  norm(X[1],2)<= 1e2
#
#X, = minimize(ls(A*x+y), [norm(x,2)<=1e2, (y-b1) in [-1e-2,1e-2]], slv)
#
#@test any(-1e-2-1e-8 .<= (X[2]-b1) .<= 1e-2+1e-8)
#@test  norm(X[1],2)<= 1e2
#
#####testing merge regularize
#X, = minimize(ls(A*x+y)+1e-3*ls(x), [norm(x,2)<=1e2, (y-b1) in [-1e-2,1e-2]], slv)
#
#@test any(-1e-2-1e-8 .<= (X[2]-b1) .<= 1e-2+1e-8)
#@test norm(X[1],2)<= 1e2
#
#X, = minimize(ls(A*x+y)+1e-3*ls(5.0*x-randn(n)), [norm(5.0*x,2)<=1e2, (y-b1) in [-1e-2,1e-2]], slv)
#
#@test any(-1e-2-1e-8 .<= (X[2]-b1) .<= 1e-2+1e-8)
#@test norm(5*X[1],2)<= 1e2

#println("testing single variable dual")
#n,m = 5,3
#x = OptVar(zeros(n))
#A = randn(m,n)
#b1 = randn(n)
#b2 = randn(m)
#d = rand(n)+2
#
#X, = minimize(ls(x-b1)+1e-10*norm(A*x-b2,1), slv)
#@test norm(X-b1) <= 1e-6
#
#r = 1e-3
#X, = minimize(ls(x-b1), norm(A*x-b2,1) <= r, slv)
#@test norm(A*X-b2,1) <= r+1e-5
#
#X, = minimize(ls(d.*x-b1)+1e-10*norm(A*x-b2,1), slv)
#@test norm(d.*X-b1) <= 1e-6
#
#r = 1e-3
#X, = minimize(ls(d.*x), norm(A*x-b2,1) <= r, slv)
#@test norm(A*X-b2,1) <= r+1e-5
#
#r = 1e8
#X, = minimize(ls(d.*x-b1), norm(A*x-b2,1) <= r, slv)
#@test norm(d.*X-b1) <= 1e-6
#
#X, = minimize(ls(x)+1e3*hingeloss(A*x,b2), slv)
#
#println("testing block variables dual")
#
#x,y = OptVar(zeros(n)),OptVar(m)
#r = 1e8
#X, = minimize(ls(d.*x-b1)+ls(5.0*y-b2), norm(A*x-y,1) <= r, slv)
#@test norm(d.*X[1]-b1) <= 1e-6
#@test norm(5.*X[2]-b2) <= 1e-6
#
#r = 1e-3
#X, = minimize(ls(d.*x-b1)+ls(5.0*y-b2), norm(A*x-y,1) <= r, slv)
#@test norm(A*X[1]-X[2],1) <= r+1e-5
#
#
