
using RegLS
using ProximalOperators
using DSP
using AcFDTD

srand(123)
X = 0.3                   # spatial sampling
env = FDTDEnv(X,SLF())

ξ = [50.;50.;10.;10.;50.;10.]; 
geo = CuboidRoom(10, 15, 12, ξ, env)
Nxyz = geo.Nx*geo.Ny*geo.Nz
Ks,Km = 1,10 
xs = [randperm(geo.Nx)[1:Ks] randperm(geo.Ny)[1:Ks] randperm(geo.Nz)[1:Ks]]'
xr = [randperm(geo.Nx)[1:Km] randperm(geo.Ny)[1:Km] randperm(geo.Nz)[1:Km]]'
indxr = sub2ind((geo.Nx,geo.Ny,geo.Nz),xr[1,:],xr[2,:],xr[3,:])
indxs = sub2ind((geo.Nx,geo.Ny,geo.Nz),xs[1,:],xs[2,:],xs[3,:])
indxr0 = collect(1:Nxyz)
deleteat!(indxr0,sort(indxr))

Nt = 600
s_true = zeros(Nt)
s_true[1] = 1.
f2 = geo.env.Fs/2*0.2#cut-off frequency of source
#filt!(s_true,digitalfilter(Bandpass(10,f2;fs = geo.env.Fs),Butterworth(5)),s_true)
s_true = repmat(s_true,1,Ks)
#[s_true[:,i] = fftfilt(s_true[:,i],randn(Nt)) for i = 1:Ks]

Qm,A,Qp = get_QmAQp(geo)
p_true = fdtd(s_true,xs,xr,Nt,geo,Qm,A,Qp)

Xs = zeros(Int64,3,Nxyz) 
Xs[1,:],Xs[2,:],Xs[3,:] = ind2sub((geo.Nx,geo.Ny,geo.Nz),collect(1:Nxyz))
#Xs = xs

L = X-> fdtd(X,Xs,xr,Nt,geo,Qm,A,Qp)
Ladj = Y-> fdtdAdj(Y,Xs,xr,Nt,geo,Qm,A,Qp)

###test adjoint
X = randn(Nt,size(Xs,2))
Y = randn(Nt,Km)
norm(vecdot(L(X),Y)-vecdot(X,Ladj(Y))) #verify adjoint operator

lambda_max = vecnorm(sqrt(sum(abs2(Ladj(p_true)),1)), Inf)  #maximum λ
lambda_max = vecnorm(Ladj(p_true), Inf)  #maximum λ
X = zeros(X)
g = NormL1(0.99*lambda_max)

solve!(L,Ladj, p_true, g, X, ZeroFPR(verbose=2,tol = 1e-10))




