using RegLS
using AbstractOperators
using RIM

srand(123)

Fs = 4000        # sampling frequency
SNR = 15
env = AcEnv(Fs)

Lx,Ly,Lz = 4.,5.,3.
T60 = 0.3
geo = CuboidRoom(Lx,Ly,Lz,T60,env)

xs = [0.5 0.5 0.5]'                    #src pos (in meters)
xr = [Lx-0.1 Ly-0.3 Lz-0.2]' #mic pos
Nh = div(Fs,4)                 #IR samples 
Nx = div(Fs,2)                 #x  samples 

h = rim(xs,xr,Nh,geo,env)[:]

x = full(sprandn(Nx, 0.05))

y0 = conv(x,h)
y = y0+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0))

H = Conv(Float64,size(x),h)

x0 = Variable(zeros(x))

lambda = 1e-2*vecnorm(H'*y,Inf)

@minimize ls(H*x0-y)+lambda*norm(x0,1)

x1 = copy(~x0)

~x0 .= 0.

@minimize ls(H*x0-y)+1e-12*norm(x0,1) with ZeroFPR(verbose = 0) 

xu = copy(~x0) #unregularized

println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

t = linspace(0,1/Fs*Nx,Nx)

using PyPlot
figure()
plot(t,x )
plot(t,x1)
plot(t,xu)

T = hcat([[zeros(i);h;zeros(Nx-1-i)] for i = 0:Nx-1]...);

~x0 .= 0.;
slv = @minimize ls(conv(x0,h)-y)+lambda*norm(x0,1); 

println("abstract operator time: $(slv.time)")

~x0 .= 0.;
slv = @minimize ls(T*x0-y)+lambda*norm(x0,1); 

println("full matrix       time: $(slv.time)")





