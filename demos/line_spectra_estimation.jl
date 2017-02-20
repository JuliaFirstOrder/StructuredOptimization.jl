
using RegLS
using ProximalOperators
using DSP

srand(12)

fs = 44e3
Nt = round(Int64,fs/35) #time samples
y0 =  zeros(Nt)

t = 0:1/fs:(Nt-1)/fs    # time axis
f = 0:fs/Nt:fs/2        # frequency axis
K = 7                   # number of sinusoids
fk = f[randperm(50)[1:K]]+f[2]*rand(K) #off axis, sinusoids frequencies
ak = randn(K)           # amplitude

for i in eachindex(fk) y0+= ak[i]*sin(2*π*fk[i].*t) end
y = y0+0.4*randn(Nt)  # sum of sinusoids corrupted by noise
Y = rfft(y)          # dft of y

zp = 10 #zeropad samples
fzp = 0:fs/(zp*Nt-1):fs
Yzp = fft([y;zeros((zp-1)*length(y))])


s = 10 #super-resolution factor
f_s = 0:fs/(s*Nt-1):fs        # super resolution frequency axis

Lifft = plan_ifft(ones(s*Nt))
Lfft  = plan_fft(ones(s*Nt))

L =    X-> (Lifft*X)[1:Nt]*sqrt(s*Nt)              #linear dft operator
Ladj = x-> Lfft*([x;zeros((s-1)*length(x))])/sqrt(s*Nt) #adjoint

x = OptVar(randn(s*Nt)+im*randn(s*Nt))
A = ifft(x)[1:Nt] 

XX = randn(s*Nt)+im*randn(s*Nt)
YY = randn(Nt) 

norm(vecdot(L(XX),YY)-vecdot(XX,Ladj(YY))) #verify adjoint operator

lambda_max = norm(Ladj(y), Inf)

X0 = zeros(Complex{Float64},s*Nt) #initial guess 
@time X,  = solve(A, y, NormL1(lambda_max*0.03), X0, ZeroFPR(tol=1e-5,verbose = 1))
Xnnz = countnz(abs(X).>=0.1*maximum(abs(X)))               #count non zero values with threshold
@time XL0,  = solve(A, y, IndBallL0(Xnnz), X, ZeroFPR(tol = 1e-5)) #solve with IndBallL0

using PyPlot
figure()
subplot(2,1,1)
plot(t,y, label = "y")
plot(t,y0, label = "ground truth")
plot(t,real(L(X)), "k", label = "recovered")
legend()
subplot(2,1,2)
plot(fzp,abs(Yzp./Nt), label = "dft zero pad.")
plot(f,abs(Y./Nt), label = "dft")
plot(fk,abs(ak)/2, "rd", label = "true amp.")
plot(f_s,abs(X)/sqrt(s*Nt), "k*", label = "LASSO")
plot(f_s,abs(XL0)/sqrt(s*Nt), "bo", label = "IndBallL0")
xlim([0;2000])
legend()




