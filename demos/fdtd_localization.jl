
using RegLS
using ProximalOperators
using DSP
using AcFDTD

srand(0)
X = 0.2                   # spatial sampling
env = FDTDEnv(X,IWB())

ξ = [100.;100.;100.;100.;100.;100.]; 
geo = CuboidRoom(2.5, 4.5, 3., ξ, env)
Nxyz = geo.Nx*geo.Ny*geo.Nz
Ks,Km = 1,4 #numbero of sources and mics
xs = [randperm(geo.Nx)[1:Ks] randperm(geo.Ny)[1:Ks] randperm(geo.Nz)[1:Ks]]'
#source position
xr = [randperm(geo.Nx)[1:Km] randperm(geo.Ny)[1:Km] randperm(geo.Nz)[1:Km]]'
#mic position
indxr = sub2ind((geo.Nx,geo.Ny,geo.Nz),xr[1,:],xr[2,:],xr[3,:])
indxs = sub2ind((geo.Nx,geo.Ny,geo.Nz),xs[1,:],xs[2,:],xs[3,:])
indxr0 = collect(1:Nxyz)
deleteat!(indxr0,sort(indxr)) #delete position where mics are

Nt = 1000
#create sound source signal
s_true = zeros(Nt)
s_true[1] = 1.
f2 = geo.env.Fs/2*0.5#cut-off frequency of source
filt!(s_true,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(5)),s_true)
s_true = repmat(s_true,1,Ks)
#[s_true[:,i] = fftfilt(s_true[:,i],randn(Nt)) for i = 1:Ks] #convolve with random noise 

Qm,A,Qp = get_QmAQp(geo)
p_true = fdtd(s_true,xs,xr,Nt,geo,Qm,A,Qp)
p_true ./= vecnorm(p_true)

Xs = zeros(Int64,3,Nxyz) 
Xs[1,:],Xs[2,:],Xs[3,:] = ind2sub((geo.Nx,geo.Ny,geo.Nz),collect(1:Nxyz))
#Xs = xs

alpha = 1e-4 #avoids too small Lip. constant
L = X-> fdtd(alpha.*X,Xs,xr,Nt,geo,Qm,A,Qp)
Ladj = Y-> fdtdAdj(Y,Xs,xr,Nt,geo,Qm,A,Qp).*alpha

###test adjoint
X = randn(Nt,size(Xs,2))
Y = randn(Nt,Km)
norm(vecdot(L(X),Y)-vecdot(X,Ladj(Y))) #verify adjoint operator

lambda_max = vecnorm(sqrt(sum(abs2(Ladj(p_true)),1)), Inf)  #maximum λ
#lambda_max = vecnorm(Ladj(p_true), Inf)  #maximum λ
X = zeros(X)
#g = NormL1(0.7*lambda_max)
g = NormL21(0.8*lambda_max)
g = IndBallL02(1,Inf)

X, = solve!(L,Ladj, p_true, g, X, ZeroFPR(verbose=2,tol = 1e-5))

using PyPlot
figure()
plot(sum(abs(X),1)')
plot(indxs-1,ones(length(indxs)).*maximum(abs(sum(abs(X),1))'), "r*")   #sources
plot(indxr-1,ones(length(indxr)).*maximum(abs(sum(abs(X),1))')/2, "bo") #mics

figure()
plot(p_true[:])
plot(L(X)[:])



