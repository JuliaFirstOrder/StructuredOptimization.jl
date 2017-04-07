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

X = Variable(N*n*m,3)
Y = Variable(N,3*n*m)
slv = ZeroFPR()
slv = FPG()
slv = slv(verbose = 1, tol = 1e-3, linesearch = false, gamma = 0.5)

@time slv = minimize(ls(reshape(X,N,3*n*m)+Y-F)+0.1*sum(norm(X),2), [rank(Y) <= 1], slv)
show(slv)

Frg = copy(~X)
Frg[Frg.!=0] = Frg[Frg.!=0]+reshape(~Y,N*m*n*3)[Frg.!=0]
Frg[repmat(sum(Frg,2).==0,1,3)] = 1 #put white in null pixels
Frg = reshape(Frg,N,3*n*m)

f = 1	
a = @sprintf("%5.5i",frames[f])
img = load("bootstrap/b$a.bmp")
Bkg = reshape((~Y)[1,:],n,m,3)
Frgf = reshape(Frg[f,:],n,m,3)
Bkgim = [ RGB(Bkg[i,ii,1],Bkg[i,ii,2],Bkg[i,ii,3]) for i = 1:n,ii =1:m]
Frgim = [ RGB(Frgf[i,ii,1],Frgf[i,ii,2],Frgf[i,ii,3]) for i = 1:n,ii =1:m]

imshow(img)
imshow(Frgim)
imshow(Bkgim)

