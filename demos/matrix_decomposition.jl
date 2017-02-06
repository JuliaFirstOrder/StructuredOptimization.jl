using RegLS
using ProximalOperators
using Images, ImageView

srand(123)
N = 20
Nframes = 2621 #number of frames
frames =  randperm(Nframes)[1:N]
n,m = 120,160 #frame size

# data from http://research.microsoft.com/en-us/um/people/jckrumm/wallflower/testimages.htm
# Bootstrapt.zip
R,G,B = zeros(Float64,n,m,N),zeros(Float64,n,m,N),zeros(Float64,n,m,N)

for f in eachindex(frames)
	a = @sprintf("%5.5i",frames[f])
	img = load("bootstrap/b$a.bmp")
	for i =1:n,ii =1:m
		R[i,ii,f] = convert(Float64,img[i,ii].r)
		G[i,ii,f] = convert(Float64,img[i,ii].g)
		B[i,ii,f] = convert(Float64,img[i,ii].b)
	end
end

F = [reshape(R,n*m,N); reshape(G,n*m,N); reshape(B,n*m,N)]'


Frg = zeros(Float64,N*n*m,3)
Bkg = zeros(Float64,N,3*n*m)
L = X-> reshape(X[1],N,3*n*m)+X[2]
Ladj = Y-> [reshape(Y,N*n*m,3),Y]

#XX =[randn(size(Frg)),randn(size(Bkg))] 
#norm(vecdot(L(XX),F)-vecdot(XX,Ladj(F))) #verify adjoint operator

g = SeparableSum([NormL21(0.1,2), IndBallRank(1)])

@time X,  = solve(L,Ladj, F, g, [Frg,Bkg], ZeroFPR(verbose = 2, tol = 1e-3))

Frg = copy(X[1])
Frg[Frg.!=0] = Frg[Frg.!=0]+reshape(X[2],N*m*n*3)[Frg.!=0]
Frg[repmat(sum(Frg,2).==0,1,3)] = 1 #put white in null pixels
Frg = reshape(Frg,N,3*n*m)

f = 1	
a = @sprintf("%5.5i",frames[f])
img = load("bootstrap/b$a.bmp")
Bkg = reshape(X[2][1,:],n,m,3)
Frgf = reshape(Frg[f,:],n,m,3)
Bkgim = [ RGB(Bkg[i,ii,1],Bkg[i,ii,2],Bkg[i,ii,3]) for i = 1:n,ii =1:m]
Frgim = [ RGB(Frgf[i,ii,1],Frgf[i,ii,2],Frgf[i,ii,3]) for i = 1:n,ii =1:m]

imshow(img)
imshow(Frgim)
imshow(Bkgim)









