verb = 1
slv = ZeroFPR(verbose = verb)

x,y,z,u = Variable(randn(4)), Variable(randn(4)), Variable(randn(2,2)), Variable(randn(4))
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

out1 = RegLS.affine(smooth_exp)[1](~x_sorted)
out2 = RegLS.affine(smooth)[1]([~x,~y])
@test norm(out1-out2) <= 1e-8

out1 = RegLS.affine(smooth_exp)[2](~x_sorted)
out2 = RegLS.affine(smooth)[2](~[y,z,u])
@test norm(out1-out2) <= 1e-8

####test single block of variable
println("testing single variable primal")
M = randn(50,50)
b = randn(50)
x = Variable(zeros(50))
P = problem(ls(M*x-b), norm(x,1)<=10)
show(P)

slv1 = solve(P,slv)
show(slv1)
@test norm(~x,1)-10 <= 1e-5
slv1 = minimize(ls(M*x-b), norm(x,1)<=10, slv)


n,m = 5,3
A = randn(m,n)
x = Variable(zeros(n))
b1,b2 = randn(m),randn(n)
minimize(ls(A*x-randn(m))+1e-2*norm(x,1), slv)

x = Variable(zeros(n))
minimize(ls(x-b2)+1e-2*norm(x,1), slv)

minimize(ls(A*x-b1), [norm(5.0*x,1) <= 1/1e-2], slv)
@test norm(~x,1) <= 1/1e-2

####test merge regularize
minimize(ls(A*x-b1)+1e7*ls(x-b2), [norm(x,1) <= 1/1e-2], slv)
@test norm(~x-b2) < 1e-5

minimize(ls(A*x-b1)+1e-7*ls(x-b2), [norm(x,1) <= 1/1e-2], slv)
@test norm(A*~x-b1) < 1e-4


#println("testing block variables primal")
n,m = 5,3
A = randn(m,n)
b1,b2 = randn(m),randn(n)
x, y = Variable(n), Variable(Complex{Float64},n)
minimize(ls(fft(x)+y-fft(b2))+1e-3*norm(x,1), slv)
cf = ls(fft(x)+y-fft(b2))+1e-3*norm(x,1)
@test norm(fft(~x)+~y-fft(b2))<1e-3

x, y = Variable(n), Variable(m)
minimize(ls(A*x+y)+1e-2*norm(x,1), norm(y-b1,1) <= 1e-2, slv)

@test (norm(~y-b1,1)-1e-2)<1e-8

minimize(ls(A*x+y), [norm(x,2)<=1e2, norm(2*y-b1,1) <= 1e-2], slv)

@test (norm(2.*~y-b1,1)-1e-2)<1e-8
@test  norm(~x,2)<= 1e2

minimize(ls(A*x+y), [norm(x,2)<=1e2, (y-b1) in [-1e-2,1e-2]], slv)

@test any(-1e-2-1e-8 .<= (~y-b1) .<= 1e-2+1e-8)
@test  norm(~x,2)<= 1e2

####testing merge regularize
n,m = 5,3
A = randn(m,n)
b1,b2 = randn(m),randn(n)
x, y = Variable(n), Variable(m)

lb,ub = -1e-5,1e-5
minimize(ls(A*x+y)+1e-3*ls(x), [norm(x,2)<=1e2, (y-b1) in [lb,ub]], slv)

@test any(lb-1e-8 .<= (~y-b1) .<= ub+1e-8)
@test norm(~x,2)<= 1e2

x, y = Variable(n), Variable(m)
minimize(ls(A*x+y)+1e-3*ls(5.0*x-randn(n)), [norm(5.0*x,2)<=1e2, (y-b1) in [lb,ub]], slv)

@test any(lb-1e-8 .<= (~y-b1) .<= ub+1e-8)
@test norm(5*~x,2)<= 1e2

#println("testing single variable dual")
#n,m = 5,3
#x = Variable(zeros(n))
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
#x,y = Variable(zeros(n)),Variable(m)
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
