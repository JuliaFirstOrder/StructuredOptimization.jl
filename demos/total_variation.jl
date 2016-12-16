
using RegLS
using ProximalOperators
using Images, ImageView
srand(123)

using TestImages
img = testimage("lena_gray")

R = zeros(size(img.data))

for i in eachindex(img.data)
	R[i] = img.data[i]
end

Nx,Ny = size(R,1),size(R,2)

dx = X-> [l == Nx ? X[l,m]-X[l-1,m] : X[l,m]-X[l+1,m] for l = 1:Nx, m = 1:Ny] 
dy = X-> [m == Ny ? X[l,m]-X[l,m-1] : X[l,m]-X[l,m+1] for l = 1:Nx, m = 1:Ny]
dxa = Y-> [l==1 ? Y[l,m] : (l==Nx-1 ? Y[l,m]-Y[l-1,m]-Y[l+1,m] : Y[l,m]-Y[l-1,m] ) 
	   for l = 1:Nx, m = 1:Ny]
dya = Y-> [m==1 ? Y[l,m] : (m==Ny-1 ? Y[l,m]-Y[l,m-1]-Y[l,m+1] : Y[l,m]-Y[l,m-1] ) 
	   for l = 1:Nx, m = 1:Ny]

L = X-> -[dx(X)[:] dy(X)[:]]
Ladj = Y-> -dxa(reshape(Y[:,1],Nx,Ny))-dya(reshape(Y[:,2],Nx,Ny))

X = randn(Nx,Ny)
Y = randn(Nx*Ny,2)

norm(vecdot(L(X),Y)-vecdot(X,Ladj(Y))) #verify adjoint operator

R_w = R+sqrt(0.006*norm(R[:],Inf))*randn(Nx,Ny)

# l-1 norm
lambda_max = vecnorm(L(R),Inf)   #this is not right...
lambda = 0.08*lambda_max
g = IndBallLinf(lambda)

### l-2 norm
#lambda_max = 100                #this is not right...
#lambda = 0.15*lambda_max
#g = IndBallL2(lambda)
#
### with group sparsity
#lambda_max = vecnorm(sqrt(sum(abs2(L(R)),1)), 1)  #this is not right...
#lambda = 0.3*lambda_max
#g = Conjugate(NormL21(lambda))

Y = 0*Y
Y, = solve(Ladj, L, R_w, g, Y, ZeroFPR(verbose = 1, tol = 1e-7))

R_r = (-(Ladj(Y)-R_w))

ImageView.view(R,xy=["y","x"])
ImageView.view(R_w,xy=["y","x"])
ImageView.view(R_r,xy=["y","x"])

