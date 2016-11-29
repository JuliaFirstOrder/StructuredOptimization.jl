using RegLS
using ProximalOperators
using Images, ImageView

srand(123)
N = 30
Nframes = 2621 #number of frames
frames =  randperm(Nframes)[1:N]
n,m = 160,120 #frame size

R,G,B = zeros(Float64,n,m,N),zeros(Float64,n,m,N),zeros(Float64,n,m,N)

for f in eachindex(frames)
	a = @sprintf("%5.5i",frames[f])
	img = data(load("bootstrap/b$a.bmp"))
	for i =1:n,ii =1:m
		R[i,ii,f] = convert(Float64,img[i,ii].r)
		G[i,ii,f] = convert(Float64,img[i,ii].g)
		B[i,ii,f] = convert(Float64,img[i,ii].b)
	end
end

RGB = [reshape(R,n*m,N); reshape(G,n*m,N); reshape(B,n*m,N)]'

Frg = zeros(Float64,N*n*m,3)
Bkg = zeros(Float64,N,3*n*m)
L = X-> reshape(X[1],N,3*n*m) +X[2]
Ladj = Y-> [reshape(Y,N*n*m,3),Y]

#norm(vecdot(L([Frg,Bkg]),RGB)-vecdot([Frg,Bkg],Ladj(RGB))) #verify adjoint operator

g = SeparableSum([NormL21(0.1,2), IndBallRank(1)])

@time X,  = solve(L,Ladj, RGB, g, [Frg,Bkg], ZeroFPR(verbose = 2))

Frg = copy(X[1])
Frg[repmat(sum(Frg,2).==0,1,3)] = 1 #put white in null pixels
Frg[Frg.<0] = abs(Frg[Frg.<0])
Frg = reshape(Frg,N,3*n*m)
Bkg = X[2][1,:]

f = 1	
a = @sprintf("%5.5i",frames[f])
img = data(load("bootstrap/b$a.bmp"))
Bkgim = shareproperties(img, reshape(Bkg,n,m,3))
Frgf = reshape(Frg[f,:],n,m,3)

Frgim = shareproperties(img, Frgf)
ImageView.view(img,xy=["y","x"])
ImageView.view(Frgim,xy=["y","x"])
ImageView.view(Bkgim,xy=["y","x"])
