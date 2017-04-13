println("\n test Variable \n")
n,m,l,i = 4,3,5,4
M = randn(n,m)
N = randn(n,l)
I = randn(n,i)
x = Variable(randn(m))
y = Variable(randn(l))
z = Variable(randn(i))
b = randn(n)
show(x)

xc = Variable(Complex{Float64},3)
show(xc)

println("\n test Affine \n")

#wrong size of L
@test_throws Exception Eye(4)*x
#wrong size of x
@test_throws Exception Eye(2)*x

A = 10*x
@test norm(A(~x)-10*~x) <= 1e-8

A = M*x
show(A)
@test norm(A(~x)-M*~x) <= 1e-8
A = M*x+M*x
@test norm(A(~x)-2*M*~x) <= 1e-8
A = M*x-M*x
@test norm(A(~x)) <= 1e-8
A = M*x+N*y
@test norm(A(~[x,y])-(M*~x+N*~y)) <= 1e-8
A = M*x-N*y
@test norm(A(~[x,y])-(M*~x-N*~y)) <= 1e-8
A = M*x-N*y+M*x
@test norm(A(~[x,y])-(2*M*~x-N*~y)) <= 1e-8
A = M*x-N*y-N*y
@test norm(A(~[x,y])-(M*~x-2*N*~y)) <= 1e-8
A = M*x+N*y+I*z
@test norm(A(~[x,y,z])-(M*~x+N*~y+I*~z)) <= 1e-8
show(A)

A = z+z
@test norm(A(~z)-(~z+~z) )<1e-9
A = z-z
@test norm( A(~z) )<1e-9
A = z+I*z
@test norm(A(~z)-(~z+I*~z) )<1e-9
A = -z+I*z
@test norm(A(~z)-(-~z+I*~z) )<1e-9
A = M*x+z
@test norm(A(~[x,z])-(M*~x+~z) )<1e-9
A = M*x-z
@test norm(A(~[x,z])-(M*~x-~z) )<1e-9


println("\n test Tilted Affine \n")

#wrong size of b
@test_throws Exception Eye(3)*x+randn(4)
@test_throws Exception Eye(3)*x+3
##wrong type of b
@test_throws Exception Eye(3)*x+(randn(4)+im)

A = 10*(z-b)
@test norm(A(~z)-(10*(~z-b)) ) < 1e-9

A = M*x+b
@test norm(A(~x)-(M*~x+b)) < 1e-9
A = M*x+N*y+b
@test norm(A(~[x,y])-(M*~x+N*~y+b)) < 1e-9
A = M*x+b+N*y
@test norm(A(~[x,y])-(M*~x+N*~y+b)) < 1e-9
A = b+M*x
@test norm(A(~x)-(M*~x+b)) < 1e-9
A = b+M*x+N*y
@test norm(A(~[x,y])-(M*~x+N*~y+b)) < 1e-9
A = M*x-b
@test norm(A(~x)-(M*~x-b)) < 1e-9
A = M*x+N*y-b
@test norm(A(~[x,y])-(M*~x+N*~y-b)) < 1e-9
A = M*x-b+N*y
@test norm(A(~[x,y])-(M*~x+N*~y-b)) < 1e-9
A = b-M*x
@test norm(A(~x)+(M*~x-b)) < 1e-9
A = b-M*x-N*y
@test norm(A(~[x,y])+(M*~x+N*~y-b)) < 1e-9
A = (M*x-b)+b
@test norm(A(~x)-(M*~x)) < 1e-9
A = b+(M*x-b)
@test norm(A(~x)-(M*~x)) < 1e-9
A = (M*x+b)-b
@test norm(A(~x)-(M*~x)) < 1e-9
A = b-(M*x+b)
@test norm(A(~x)+(M*~x)) < 1e-9

show(A)
Z = Variable(randn(4,4))

println("\n testing composition affine \n")

A = DFT(n)*(M*x)
@test norm(A(~x)-fft(M*~x)) < 1e-9

A = DFT(n)*(M*x-b)
@test norm(A(~x)-fft(M*~x-b)) < 1e-9

A = pi*(M*x)
@test norm(A(~x)-pi*(M*~x)) < 1e-9

A = 4*(M*x-b)
@test norm(A(~x)-4*(M*~x-b)) < 1e-9


println("\n testing function constructors \n")
show(zeros(z)-b)
show(eye(z)-b)
show(reshape(z,2,2)-reshape(b,2,2))
show(randn(n).*z-b)
show(randn(n,n)*z-b)
show(z[1:2]-b[1:2])
show(fft(z))
show(ifft(z))
show(dct(z))
show(idct(z))
show(finitediff(z))
show(finitediff(Z,2))
show(variation(Z))
show(conv(z,randn(10)))
show(xcorr(z,randn(10)))

println("\n testing function composition \n")
A = zeros(M*x-b)
@test norm(A(~x)) <1e-9

A = eye(M*x)
@test norm(A(~x)-(M*~x)) <1e-9

A = reshape(M*x,2,2)
@test norm(A(~x)-reshape(M*~x,2,2)) <1e-9

A = reshape(M*x-b,2,2)
@test norm(A(~x)-reshape(M*~x-b,2,2)) <1e-9

A = (2*ones(n)).*(M*x-b)
@test norm(A(~x)-2*(M*~x-b)) <1e-9

A = (M*x)[1:2]
@test norm(A(~x)-(M*~x)[1:2]) <1e-9

A = (4*ones(5,4))*(M*x-b)
@test norm(A(~x)-(4*ones(5,4))*(M*~x-b)) <1e-9

A =  fft(M*x-b)
@test norm(A(~x)-fft(M*~x-b)) <1e-9

A =  ifft(M*x-b)
@test norm(A(~x)-ifft(M*~x-b)) <1e-9

A =  dct(M*x-b)
@test norm(A(~x)-dct(M*~x-b)) <1e-9

A =  finitediff(M*x-b)
@test norm(A(~x)-FiniteDiff((n,))*(M*~x-b)) <1e-9

X = Variable(randn(3,6))
B = randn(4,6)
A =  finitediff(M*X-B,2)
@test norm(A(~X)-FiniteDiff((4,6),2)*(M*~X-B)) <1e-9

A =  variation(M*X-B)
@test norm(A(~X)-Variation(4,6)*(M*~X-B)) <1e-9

A =  conv(M*x-b,b)
@test norm(A(~x)-conv(M*~x-b,b)) <1e-9

A =  xcorr(M*x-b,b)
@test norm(A(~x)-xcorr(M*~x-b,b)) <1e-9


println("\n testing sort_and_expand Affine \n")

##test sort_and_expand

x,y,z = Variable(10), Variable(2,5), Variable(10)
X = [randn(10), randn(2,5), randn(10)]
Y = randn(10)
A = reshape(y,10)

B = RegLS.sort_and_expand([x,y,z],A)

test1,test2 = RegLS.test_FwAdj(operator(B), X, Y)
@test test1 < 1e-8
@test test2 < 1e-8
test3 = RegLS.test_Op(operator(B), X, Y)
@test test3 < 1e-8

A = reshape(y,10)+2.5*x+randn(10)

B = RegLS.sort_and_expand([x,y,z],A)

test1,test2 = RegLS.test_FwAdj(operator(B), X, Y)
@test test1 < 1e-8
@test test2 < 1e-8
test3 = RegLS.test_Op(operator(B), X, Y)
@test test3 < 1e-8

A = reshape(y,10)+2.5*x+randn(10,10)*z

B = RegLS.sort_and_expand([x,y,z],A)

test1,test2 = RegLS.test_FwAdj(operator(B), X, Y)
@test test1 < 1e-8
@test test2 < 1e-8
test3 = RegLS.test_Op(operator(B), X, Y)
@test test3 < 1e-8












