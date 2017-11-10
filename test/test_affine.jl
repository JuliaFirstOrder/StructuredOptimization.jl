@printf("\nTesting variables\n")

n, m, k = 3, 4, 5
x1t = Variable(Float32, n)
x1 = Variable(n)
x1i = Variable(randn(n))
x2 = Variable(n, m)
x2t = Variable(Float32, n, m)
x2i = Variable(randn(n, m))
x3 = Variable(n, m, k)
x3t = Variable(Float32, n, m, k)
xx = randn(n,m,k)
x3i = Variable(xx)
@test xx == (~x3i)

@printf("\nTesting linear expressions\n")

#### * ####
n, m1, m2, k = 3, 4, 5, 6
x1 = Variable(randn(m1))
x2 = Variable(randn(m2))
A1 = randn(n, m1)
A2 = randn(n, m2)
b  = randn(n)
b0  = pi
opA1 = MatrixOp(A1)
opA2 = MatrixOp(A2)
# multiply with Variable
ex1 = opA1*x1 
@test variables(ex1) == (x1,)
@test operator(ex1)*(~x1) == A1*(~x1)
B1 = randn(k, n)
opB1 = MatrixOp(B1)
# multiply with Expression
ex2 = opB1*ex1
@test variables(ex2) == (x1,)
@test norm(operator(ex2)*(~x1) - B1*A1*(~x1)) < 1e-12
# multiply with Expression with multiple variables
ex3 = opB1*(opA1*x1+opA2*x2)
@test variables(ex3) == (x1,x2)
@test norm(operator(ex3)*(~x1,~x2) - B1*(A1*(~x1)+A2*(~x2))) < 1e-12
# multiply with displacemented Array Expression with multiple variables
ex3 = opB1*(opA1*x1+opA2*x2+b)
@test variables(ex3) == (x1,x2)
@test norm(displacement(ex3) - B1*b) < 1e-12
# multiply with displacemented scalar Expression with multiple variables
ex3 = opB1*(opA1*x1+opA2*x2+b0)
@test variables(ex3) == (x1,x2)
@test norm(displacement(ex3) - B1*(ones(size(B1,2))*b0)) < 1e-12
@test_throws ArgumentError MatrixOp(randn(n,m1+1))*x1
@test_throws ArgumentError MatrixOp(randn(n,m1))*Variable(Complex{Float64},m1)

n, m1, m2, k = 3, 4, 5, 6
A1 = randn(n, m1)
A2 = randn(n, m2)
b1 = randn(n,n)
b2 = randn(n,n)
opA1 = MatrixOp(A1,n)
opA2 = MatrixOp(A2,n)
x1, x2  = Variable(randn(m1,n)), Variable(randn(m2,n))
# multiply Expressions 
ex = (opA1*x1)*(opA2*x2)
@test variables(ex) == (x1,x2)
@test norm(operator(ex)*(~x1,~x2) - (A1*(~x1))*(A2*(~x2))) < 1e-12
ex = (opA1*x1-b1)*(opA2*x2)
@test variables(ex) == (x1,x2)
@test norm(operator(ex)*(~x1,~x2) - (A1*(~x1)-b1)*(A2*(~x2))) < 1e-12
ex = (opA1*x1)*(opA2*x2+b2)
@test variables(ex) == (x1,x2)
@test norm(operator(ex)*(~x1,~x2) - (A1*(~x1))*(A2*(~x2)+b2)) < 1e-12
ex = (opA1*x1+b1)*(opA2*x2+b2)
@test variables(ex) == (x1,x2)
@test norm(operator(ex)*(~x1,~x2)+displacement(ex) - (A1*(~x1)+b1)*(A2*(~x2)+b2)) < 1e-12

##### reshape ####
m,n = 8,10
A = randn(n,m)
x = Variable(randn(m))
b = randn(n)
B = reshape(b,2,5)

ex = reshape(x,4,2)
@test vecnorm(operator(ex)*~x - reshape(~x,4,2)) < 1e-12

ex2 = reshape(A*x,2,5)
@test vecnorm(operator(ex2)*~x - reshape(A*~x,2,5)) < 1e-12

ex3 = reshape(A*x,2,5)+B
@test vecnorm(operator(ex2)*~x+displacement(ex3)- reshape(A*~x,2,5)-B) < 1e-12

ex4 = reshape(A*x-b,2,5)
@test vecnorm(operator(ex4)*~x+displacement(ex4)- reshape(A*~x-b,2,5)) < 1e-12

#### + ####

# sum same variable
n, m = 3,4
x = Variable(randn(m))
A = randn(n, m)
B = randn(n, m)
opA = MatrixOp(A)
opB = MatrixOp(B)
ex1 = opA*x+opB*x 
@test variables(ex1) == (x,)
@test norm(operator(ex1)*~x - (opA+opB)*~x) < 1e-12

# sum different variables
n, m1, m2 = 3, 4, 5
xa = Variable(randn(m1))
xb = Variable(randn(m2))
xc = Variable(randn(n))
xd = Variable(randn(n))
A = randn(n, m1)
A2 = randn(n, m1)
B = randn(n, m2)
b = randn(n)
b0 = pi
opA = MatrixOp(A)
opA2 = MatrixOp(A2)
opB = MatrixOp(B)
opI = Eye(n)

# (+) sum different variables no HCAT
ex1 = opA*xa+opB*xb 
@test variables(ex1) == (xa,xb)
@test norm(operator(ex1)*(~xa,~xb) - hcat(opA,opB)*(~xa,~xb)) <1e-12

# (+) sum of same variables
ex2 = opA2*xa 
exs1 = ex1+ex2
exs2 = ex2+ex1
@test variables(exs1) == (xa,xb)
@test norm(operator(exs1)*(~xa,~xb) - hcat(opA+opA2,opB)*(~xa,~xb)) <1e-12
@test variables(exs2) == (xa,xb)
@test norm(operator(exs2)*(~xa,~xb) - hcat(opA+opA2,opB)*(~xa,~xb)) <1e-12

# (+) sum of different variables with HCAT
exs3 = exs1+exs2
@test variables(exs3) == (xa,xb)
@test norm(operator(exs3)*(~xa,~xb) - hcat(2*(opA+opA2),2*opB)*(~xa,~xb)) <1e-12
# (+) sum of different variables with HCAT
exs4 = exs1+(xc+xd)
@test variables(exs4) == (xa,xb,xc,xd)
@test norm(operator(exs4)*(~xa,~xb,~xc,~xd) - hcat(opA+opA2,opB,opI,opI)*(~xa,~xb,~xc,~xd)) <1e-12

# (+) sum Array
ex1 = xd+b
@test norm(displacement(ex1) - b) == 0.

# (+) sum scalar
ex2 = opB*xb+b0
@test (displacement(ex2) - b0) == 0.

# (+) sum displacemented expressions
ex3 = ex1+ex2
@test norm(displacement(ex3) - (b+b0)) == 0.

# (-) sum different variables no HCAT
ex1 = opA*xa-opB*xb 
@test variables(ex1) == (xa,xb)
@test norm(operator(ex1)*(~xa,~xb) - hcat(opA,-opB)*(~xa,~xb)) <1e-12

# (-) sum of same variables
ex2 = opA2*xa 
exs1 = ex1-ex2
exs2 = ex2-ex1
@test variables(exs1) == (xa,xb)
@test norm(operator(exs1)*(~xa,~xb) - hcat(opA-opA2,-opB)*(~xa,~xb)) <1e-12
@test variables(exs2) == (xa,xb)
@test norm(operator(exs2)*(~xa,~xb) - hcat(-opA+opA2,+opB)*(~xa,~xb)) <1e-12

# (-) sum of same variables with HCAT
exs3 = exs1-exs2
@test variables(exs3) == (xa,xb)
@test norm(operator(exs3)*(~xa,~xb) - hcat(2*(opA-opA2),-2*opB)*(~xa,~xb)) < 1e-12
# (-) sum of different variables with HCAT
exs4 = exs1-(xc-xd)
@test variables(exs4) == (xa,xb,xc,xd)
@test norm(operator(exs4)*(~xa,~xb,~xc,~xd) - hcat(opA-opA2,-opB,-opI,opI)*(~xa,~xb,~xc,~xd)) <1e-12

# (-) sum Array
ex1 = xd-b
@test norm(displacement(ex1) + b) == 0.

# (-) sum scalar
ex2 = opB*xb-b0
@test (displacement(ex2) + b0) == 0.

# (+) sum displacemented expressions
ex3 = ex1-ex2
@test norm(displacement(ex3) - (-b+b0)) == 0.

@test_throws ArgumentError MatrixOp(randn(10,20))*Variable(20)+randn(11)
@test_throws ArgumentError MatrixOp(randn(10,20))*Variable(20)+(3+im)

@printf("\nTesting AbstractOperators binding\n")

# MatrixOp
n,m = 3,4
A = randn(n,m)
op = MatrixOp(A)
x = Variable(randn(m))
ex = A*x
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# LMatrixOp
n,m = 3,4
b = randn(m)
op = LMatrixOp(Float64,(n,m),b)
X = Variable(randn(n,m))
ex = X*b
@test norm(operator(ex)*(~X)-op*(~X)) <1e-12

# DiagOp
n = 3
d = randn(n)
op = DiagOp(Float64,(n,),d)
x = Variable(randn(n))
ex = d.*x
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

## GetIndex
n = 5
op = GetIndex(Float64,(n,),(1:2,))
x = Variable(randn(n))
ex = x[1:2]
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# DFT
n = 5
op = DFT(Float64,(n,))
x = Variable(randn(n))
ex = fft(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# IDFT
n = 5
op = IDFT(Float64,(n,))
x = Variable(randn(n))
ex = ifft(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# RDFT
n = 5
op = RDFT(Float64,(n,))
x = Variable(randn(n))
ex = rfft(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# IRDFT
n = 5
op = IRDFT(Complex{Float64},(n,),8)
x = Variable(randn(n)+im*randn(n))
ex = irfft(x,8)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# DCT
n = 5
op = DCT(Float64,(n,))
x = Variable(randn(n))
ex = dct(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# IDCT
n = 5
op = IDCT(Float64,(n,))
x = Variable(randn(n))
ex = idct(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Conv
n = 5
h = randn(n)
op = Conv(Float64,(n,),h)
x = Variable(randn(n))
ex = conv(x,h)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Xcorr
n = 5
h = randn(n)
op = Xcorr(Float64,(n,),h)
x = Variable(randn(n))
ex = xcorr(x,h)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Filt
n = 5
h = randn(n)
op = Filt(Float64,(n,),h)
x = Variable(randn(n))
ex = filt(x,h)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# MIMOFilt
m,n = 10,2
b = [[1.;0.;1.;0.;0.],[1.;0.;1.;0.;0.]]
a = [[1.;1.;1.],[2.;2.;2.]]
op = MIMOFilt(Float64,(m,n),b,a)
x = Variable(randn(m,n))
ex = mimofilt(x,b,a)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# FiniteDiff
n,m = 5,7
op = FiniteDiff(Float64,(n,m))
x = Variable(randn(n,m))
ex = finitediff(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

op = FiniteDiff(Float64,(n,m),2)
x = Variable(randn(n,m))
ex = finitediff(x,2)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Variation
n,m = 5,7
op = Variation(Float64,(n,m))
x = Variable(randn(n,m))
ex = variation(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# ZeroPad
n = 5
op = ZeroPad(Float64,(n,),10)
x = Variable(randn(n))
ex = zeropad(x,10)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Sigmoid
n = 5
op = Sigmoid(Float64,(n,),10)
x = Variable(randn(n))
ex = sigmoid(x,10)
ex = σ(x,10)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12
