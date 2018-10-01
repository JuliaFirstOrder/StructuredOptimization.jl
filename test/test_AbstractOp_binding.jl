println("\nTesting AbstractOperators binding\n")

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

# DiagOp with Scalar
n = 3
d = randn(n)
op = DiagOp(Float64,(n,),d)
x = Variable(randn(n))
ex = d.*x
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12
ex = x.*d
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Scale
n = 3
d = 5
x = Variable(randn(n))
ex = d*x
@test norm(operator(ex)*(~x)-5*(~x)) <1e-12
ex = x*d
@test norm(operator(ex)*(~x)-5*(~x)) <1e-12

## GetIndex
n = 5
op = GetIndex(Float64,(n,),(1:2,))
x = Variable(randn(n))
ex = x[1:2]
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# DFT
n = 5
op = AbstractOperators.DFT(Float64,(n,))
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

# Pow
n = 5
op = Pow(Float64,(n,),2)
x = Variable(randn(n))
ex = pow(x,2)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Exp
n = 5
op = Exp(Float64,(n,))
x = Variable(randn(n))
ex = exp(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Cos
n = 5
op = Cos(Float64,(n,))
x = Variable(randn(n))
ex = cos(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Sin
n = 5
op = Sin(Float64,(n,))
x = Variable(randn(n))
ex = sin(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Atan
n = 5
op = Atan(Float64,(n,))
x = Variable(randn(n))
ex = atan(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12

# Tanh
n = 5
op = Tanh(Float64,(n,))
x = Variable(randn(n))
ex = tanh(x)
@test norm(operator(ex)*(~x)-op*(~x)) <1e-12
