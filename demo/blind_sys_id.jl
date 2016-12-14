
using RegLS
using ProximalOperators
using RIM
using DSP

srand(123)
Fs = 8e3
Nx = 4000    #length of input singal
Nh = 300    #length of IR

env = AcEnv(Fs)

Km = 3 #num of mics
Lx = [5.;6.;3.]
geo = CuboidRoom(Lx[1],Lx[2],Lx[3],0.1,env);
xs = rand(3).*Lx
xr = rand(3,Km).*Lx

h_true = rim(xs,xr,Nh,geo,env; N = [4;4;4], Tw = 20);
h_true = randn(size(h_true))

y = zeros(Nx,Km)
x = [randn(round(Int64,Nx/2));zeros(round(Int64,Nx/2))]
x = randn(Nx)
#x = convert(Array{Float64},bitrand(Nx))

for i = 1:Km 
	y[:,i] = fftfilt([h_true[:,i];zeros(size(y,1)-Nh)],x)
end

X = randn(Nh*Km)
Y = randn(Nx,round(Int64,Km*(Km-1)/2))

fxcross = (x,y) -> flipdim(fftfilt(x,flipdim(y,1)),1)[1:Nh]
function L2(X::Array{Float64},Nh::Int64,Km::Int64,y::Array{Float64})
	X = reshape(X,Nh,Km)
	out = zeros(Float64,Nx,round(Int64,Km*(Km-1)/2))
	counter = 1
	for i = 1:Km-1,j = i+1:Km
		out[:,counter] = (fftfilt([X[:,j];zeros(size(y,1)-Nh)],y[:,i])
		                 -fftfilt([X[:,i];zeros(size(y,1)-Nh)],y[:,j])  )
		counter += 1 
	end
	return out
end

function L2adj(Y::Array{Float64},Nh::Int64,Km::Int64,y::Array{Float64})
	out = zeros(Float64,Nh,Km)
	
	counter = 1
	for i = 1:Km-1,j = i+1:Km
		out[:,j] += fxcross(y[:,i],Y[:,counter])
		out[:,i] -= fxcross(y[:,j],Y[:,counter])
		counter += 1 
	end
	return out[:]
end

L = X-> L2(X,Nh,Km,y)
Ladj = Y-> L2adj(Y,Nh,Km,y)

norm(vecdot(L(X),Y)-vecdot(X,Ladj(Y))) #verify adjoint operator

h = zeros(Nh,Km)
h = h[:]

# ||x|| = 1
g = IndSphereSqrL2(1.)  

# x[ind] = 1
#indx = [round(Int64,Nh/2)]
#indx0 = collect(1:Km*Nh)
#deleteat!(indx0,indx)
#g = SlicedSeparableSum([IndPoint(1.)=>indx,IndFree()=>indx0])
h, = solve(L,Ladj, zeros(Y), g, h, ZeroFPR( tol = 1e-6))

h = reshape(h,Nh,Km)

using PyPlot

h =  circshift(h,0)
for i = 1:Km
subplot(2,Km,i)
plot(h[:,i]/norm(h[:,i]))
plot(h_true[:,i]/norm(h_true[:,i]))
subplot(2,Km,Km+i)
plot(10*log10(abs(fft(h[:,i]/norm(h[:,i])))))
plot(10*log10(abs(fft(h_true[:,i]/norm(h_true[:,i])))))
xlim([0;round(Int64,Nh/2)])
end

NPM = 10*log10( vecnorm(h_true-vecdot(h,h_true)/vecdot(h,h)*h   )^2/vecdot(h_true,h_true) )




