#test single block of variable

x = OptVar(zeros(5))
minimize(ls(fft(x)-fft(randn(5)))+1e-2*norm(x,1) )

x = OptVar(zeros(5))
minimize(ls(x-randn(5))+1e-2*norm(x,1) )

x = OptVar(zeros(5))
minimize(ls(fft(x)-fft(randn(5))), norm(x,1) <= 1/1e-2 )

#test merge regularize
x = OptVar(zeros(5))
minimize(ls(fft(x)-fft(randn(5)))+1e-7*ls(x-randn(5)), norm(x,1) <= 1/1e-2 )

##test 2 blocks of variables
x, y = OptVar(5), OptVar(Complex{Float64},5)
minimize(ls(fft(x)+y-fft(randn(5)))+1e-2*norm(x,1) )

x, y = OptVar(5), OptVar(5)
b = randn(5) 
X, = minimize(ls(dct(x)+eye(y))+1e-2*norm(x,1), norm(y-b,1) <= 1e-2  )

@test (norm(X[2]-b,1)-1e-2)<1e-8

x, y = OptVar(5), OptVar(5)
b = randn(5) 
X, = minimize(ls(dct(x)+eye(y)), norm(x,2)<=1e2, norm(diagop(y,2)-b,1) <= 1e-2  )
X, = minimize(ls(dct(x)+eye(y)), [norm(x,2)<=1e2, norm(diagop(y,2)-b,1) <= 1e-2] )

@test (norm(2.*X[2]-b,1)-1e-2)<1e-8
@test  norm(X[1],2)<= 1e2

x, y = OptVar(5), OptVar(5)
b = randn(5) 
X, = minimize(ls(dct(x)+eye(y)), norm(x,2)<=1e2, (y-b) in [-1e-2,1e-2]  )

@test any(-1e-2-1e-8 .<= (X[2]-b) .<= 1e-2+1e-8)
@test  norm(X[1],2)<= 1e2

#test merge regularize
x, y = OptVar(5), OptVar(5)
b = randn(5) 
X, = minimize(ls(dct(x)+eye(y))+1e-3*ls(x) , norm(x,2)<=1e2, (y-b) in [-1e-2,1e-2]  )

@test any(-1e-2-1e-8 .<= (X[2]-b) .<= 1e-2+1e-8)
@test  norm(X[1],2)<= 1e2


		



