using RegLS
using AbstractOperators
using DSP

srand(17)

fs = 8e3
Nt = round(Int64,fs/35) #time samples
y0 = zeros(Nt)
SNR = 15

s = 10 #super-resolution factor
f_s = 0:fs/(s*Nt-1):fs        # super resolution frequency axis
t = 0:1/fs:(Nt-1)/fs    # time axis
f = 0:fs/Nt:fs/2        # frequency axis
K = 14                  # number of sinusoids
fk = f_s[randperm(500)[1:K]] #sinusoids frequencies
ak = 0.1*randn(K)+0.7           # amplitude

for i in eachindex(fk) y0 .+= ak[i].*sin.(2*π*fk[i].*t) end
y = y0+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0))
xdft = rfft(y)          # dft of y

xzp = fft([y;zeros((s-1)*length(y))])


F = (DFT(s*Nt)')[1:Nt] 
lambda_max = norm(F'*y, Inf)

x = Variable(zeros(Complex{Float64},s*Nt))
lambda = 0.02*lambda_max

@minimize ls(F*x-y)+lambda*norm(x,1) with ZeroFPR(tol = 1e-4)

x1 = copy(~x)

#~x .= x1.*(abs.(x1) .> 0.07)

@minimize ls(F*x-y) st norm(x,0) <= K*2 with ZeroFPR(tol = 1e-4)

x0 = copy(~x)

Nf = 573 #just to cut off stuff not plotted 

using PyPlot
figure()
subplot(2,1,1)
plot(t,y, label = "y")
plot(t,y0, label = "ground truth")
plot(t,F*x1, "k", label = "recovered")
legend()
subplot(2,1,2)
plot(f_s[1:Nf],abs.(xzp./Nt )[1:Nf] , label = "dft zero pad.")
plot(f,        abs.(xdft./Nt)       , label = "dft")
plot(fk,       abs.(ak)/2     , "rd", label = "true amp.")
plot(f_s[1:Nf],abs.(x1)[1:Nf] , "k*", label = "LASSO")
#plot(f_s[1:Nf],abs.(x1t)[1:Nf] , "mo", label = "IC")
plot(f_s[1:Nf],abs.(x0)[1:Nf] , "go", label = "IndBallL0")
xlim([0;2000])
legend()

