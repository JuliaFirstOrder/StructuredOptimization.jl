##test single block of variable
#x = OptVar(zeros(5))
#minimize(ls(fft(x)-fft(randn(5)))+1e-2*norm(x,1) )
#
#n,m = 3,5
#A = randn(m,n)
#x = OptVar(zeros(n))
#minimize(ls(A*x-randn(m))+1e-2*norm(x,1) )
#
#x = OptVar(zeros(n))
#minimize(ls(x-randn(n))+1e-2*norm(x,1) )
#
#x = OptVar(zeros(n))
#minimize(ls(A*x-randn(m)), norm(x,1) <= 1/1e-2 )
#
##test merge regularize
#x = OptVar(zeros(n))
#minimize(ls(A*x-randn(m))+1e-7*ls(x-randn(n)), norm(x,1) <= 1/1e-2 )
#
###test 2 blocks of variables
#x, y = OptVar(n), OptVar(Complex{Float64},n)
#minimize(ls(fft(x)+y-fft(randn(n)))+1e-2*norm(x,1) )
#
#x, y = OptVar(n), OptVar(m)
#b = randn(m) 
#X, = minimize(ls(A*x+y)+1e-2*norm(x,1), norm(y-b,1) <= 1e-2  )
#
#@test (norm(X[2]-b,1)-1e-2)<1e-8
#
#X, = minimize(ls(A*x+y), norm(x,2)<=1e2, norm(2*y-b,1) <= 1e-2  )
#X, = minimize(ls(A*x+y), [norm(x,2)<=1e2, norm(2*y-b,1) <= 1e-2] )
#
#@test (norm(2.*X[2]-b,1)-1e-2)<1e-8
#@test  norm(X[1],2)<= 1e2
#
#X, = minimize(ls(A*x+y), norm(x,2)<=1e2, (y-b) in [-1e-2,1e-2]  )
#
#@test any(-1e-2-1e-8 .<= (X[2]-b) .<= 1e-2+1e-8)
#@test  norm(X[1],2)<= 1e2
#
##test merge regularize
#X, = minimize(ls(A*x+y)+1e-3*ls(x), norm(x,2)<=1e2, (y-b) in [-1e-2,1e-2]  )
#
#@test any(-1e-2-1e-8 .<= (X[2]-b) .<= 1e-2+1e-8)
#@test  norm(X[1],2)<= 1e2


##test dual single variable
srand(0)
n,m = 5,3
x = OptVar(zeros(n))
A = randn(m,n)
b1 = randn(n)
b2 = randn(m)

X, = minimize(ls(x-b1)+1e-10*norm(A*x,1)  )
@test norm(X-b1) <= 1e-8

X, = minimize(ls(x-b1)+1e4*norm(A*x,1)  )
@test norm(A*X,1) <= 1e-4

r = 1e-5
X, = minimize(ls(x), norm(A*x-b2,1) <= r  )
@test norm(A*X-b2,1) <= r+1e-7

X, = minimize(ls(5.0*x-b1)+1e-10*norm(A*x,1)  )
@test norm(5.0*X-b1) <= 1e-8

