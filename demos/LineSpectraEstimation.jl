using RegLS
using AbstractOperators
using DSP

srand(17)
save_stuff = false

fs = 16e3
Nt = 2^8 #time samples
y0 = zeros(Nt)
SNR = 10

s = 6 #super-resolution factor
f_s  = linspace(0,  fs,       s*Nt+1)[1:end-1]     # super resolution frequency axis
f_s2 = linspace(0,fs/2,div(s*Nt,2)+1)              # super resolution frequency axis (up to Nyquist)
t  = 0:1/fs:(Nt-1)/fs                              # time axis
f  = linspace(0,fs,Nt+1)[1:end-1]                  # frequency axis
f2 = linspace(0,fs/2,div(Nt,2)+1)                  # frequency axis
K = 14                                             # number of sinusoids
fk = f_s2[randperm(div(s*Nt,2)+1)[1:K]]            #sinusoids frequencies
ak = 0.1*randn(K)+0.7           # amplitude

for i in eachindex(fk) y0 .+= ak[i].*sin.(2*π*fk[i].*t) end
y = y0.+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0))
xdft = fft(y)          # dft of y

xzp = rfft([y;zeros((s-1)*length(y))])

F = (s*Nt*IRDFT((div(s*Nt,2)+1,),s*Nt))[1:Nt] 
lambda_max = norm(F'*y, Inf)

x = Variable(zeros(Complex{Float64},div(s*Nt,2)+1))
lambda = 0.06*lambda_max

@minimize ls(F*x-y)+lambda*norm(x,1)  
x1 = copy(~x)

@minimize ls(F*x-y) st norm(x,0) <= K 
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
plot(f_s2,abs.(xzp./Nt ), label = "dft zero pad.")
plot(f,        abs.(xdft./Nt)       , label = "dft")
plot(fk,       abs.(ak)/2     , "r*", label = "true amp.")
plot(f_s2,abs.(x1), "k*", label = "LASSO")
plot(f_s2,abs.(x0), "go", label = "IndBallL0")
xlim([0;fs/2])
legend()
