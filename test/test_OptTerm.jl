verb = 0
slv = ZeroFPR(verbose = verb)

###test single block of variable
println("testing single variable primal")
x = OptVar(zeros(5))
minimize(ls(fft(x)-fft(randn(5)))+1e-2*norm(x,1), slv)

n,m = 3,5
A = randn(m,n)
x = OptVar(zeros(n))
minimize(ls(A*x-randn(m))+1e-2*norm(x,1), slv)

x = OptVar(zeros(n))
minimize(ls(x-randn(n))+1e-2*norm(x,1), slv)

x = OptVar(zeros(n))
X, = minimize(ls(A*x-randn(m)), [norm(x,1) <= 1/1e-2], slv)
@test norm(X,1) <= 1/1e-2

#test merge regularize
x = OptVar(zeros(n))
minimize(ls(A*x-randn(m))+1e-7*ls(x-randn(n)), [norm(x,1) <= 1/1e-2], slv)

println("testing block variables primal")
x, y = OptVar(n), OptVar(Complex{Float64},n)
minimize(ls(fft(x)+y-fft(randn(n)))+1e-2*norm(x,1), slv)

x, y = OptVar(n), OptVar(m)
b = randn(m) 
X, = minimize(ls(A*x+y)+1e-2*norm(x,1), norm(y-b,1) <= 1e-2, slv)

@test (norm(X[2]-b,1)-1e-2)<1e-8

X, = minimize(ls(A*x+y), [norm(x,2)<=1e2, norm(2*y-b,1) <= 1e-2], slv)

@test (norm(2.*X[2]-b,1)-1e-2)<1e-8
@test  norm(X[1],2)<= 1e2

X, = minimize(ls(A*x+y), [norm(x,2)<=1e2, (y-b) in [-1e-2,1e-2]], slv)

@test any(-1e-2-1e-8 .<= (X[2]-b) .<= 1e-2+1e-8)
@test  norm(X[1],2)<= 1e2

# testing merge regularize
X, = minimize(ls(A*x+y)+1e-3*ls(x), [norm(x,2)<=1e2, (y-b) in [-1e-2,1e-2]], slv)

@test any(-1e-2-1e-8 .<= (X[2]-b) .<= 1e-2+1e-8)
@test  norm(X[1],2)<= 1e2


println("testing single variable dual")
srand(23)
n,m = 5,3
x = OptVar(zeros(n))
A = randn(m,n)
b1 = randn(n)
b2 = randn(m)
d = rand(n)+2

X, = minimize(ls(x-b1)+1e-10*norm(A*x-b2,1), slv)
@test norm(X-b1) <= 1e-6

r = 1e-3
X, = minimize(ls(x-b1), norm(A*x-b2,1) <= r, slv)
@test norm(A*X-b2,1) <= r+1e-5

X, = minimize(ls(diagop(x,d)-b1)+1e-10*norm(A*x-b2,1), slv)
@test norm(d.*X-b1) <= 1e-6

r = 1e-3
X, = minimize(ls(diagop(x,d)), norm(A*x-b2,1) <= r, slv)
@test norm(A*X-b2,1) <= r+1e-5

#r = 1e8
#X, = minimize(ls(diagop(x,d)-b1), norm(A*x-b2,1) <= r, slv)
#@test norm(d.*X-b1) <= 1e-6


