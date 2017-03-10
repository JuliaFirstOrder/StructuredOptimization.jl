
using RegLS
using ProximalOperators

Nx = 1000    #length of input singal
Nh = 700     #length of IR
Ny = Nx+Nh-1 #length of output

srand(123)
x_true = abs(full(sprandn(Nx,1,0.05)))[:]   #spikes
x = OptVar(Nx)
h = exp(-(0:Nh-1)./(0.15*Nh))          #exponential IR


SNR = 20                 #signal to noise ratio
A   = conv(h,x)          #convolution operator
y_clean = A*x_true
noise = sqrt(mean(y_clean.^2)*10^(-SNR/10))*randn(Ny)
y = y_clean+noise

H = zeros(Ny,Nx)       #Toeplitz matrix for linear convolution
for i = 1:Nx
	H[i:(i-1)+Nh,i] = h
end

# test operators are the same
norm(H*x_true-A*x_true)

x_LS = H\y #least suqares solution

##non-negative deconvolution with H
@time x_NN,  = solve(H, -y, IndBox(0., Inf))
##Matrix free non-negative deconvolution
@time x_NN_MF, = minimize(ls(conv(h,x)-y), x >= 0. )

using PyPlot
figure()
subplot(2,1,1)
plot(y, label = "y")
xlim([1,length(y)])
legend()
subplot(2,1,2)
plot(x_true, label = "ground truth")
plot(x_NN, label = "NND")
#plot(x_LS, label = "LS")
xlim([1,length(y)])
legend()

