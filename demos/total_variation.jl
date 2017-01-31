
using RegLS
using ProximalOperators
using Images
using ImageView
srand(123)

using TestImages
img = testimage("lena_gray")

R = zeros(size(img.data))

for i in eachindex(img.data)
	R[i] = img.data[i]
end

Nx,Ny = size(R,1),size(R,2)
Dx = spdiagm((ones(Nx),-ones(Nx-1)),(0,1),Nx,Nx)
Dx[end,end] = 1
Dy = spdiagm((ones(Ny),-ones(Ny-1)),(0,1),Ny,Ny)
Dy[end,end-1] = -1
DDx = kron(speye(Ny),Dy)
DDy = kron(Dx,speye(Ny))

L = X-> [DDx*X[:] DDy*X[:]]
Ladj = Y-> reshape((DDx'*Y[:,1]+DDy'*Y[:,2]),Nx,Ny)

X = randn(Nx,Ny)
Y = randn(Nx*Ny,2)

norm(vecdot(L(X),Y)-vecdot(X,Ladj(Y))) #verify adjoint operator

R_w = R+sqrt(0.006*norm(R[:],Inf))*randn(Nx,Ny)

# l-1 norm
lambda_max = vecnorm(L(R),Inf)   #this is not right...
lambda = 0.08*lambda_max
g = NormL1(lambda)

### l-2 norm
#lambda_max = 100                #this is not right...
#lambda = 0.1*lambda_max
#g = NormL2(lambda)
##
#### with group sparsity
#lambda_max = vecnorm(sqrt(sum(abs2(L(R)),1)), 1)  #this is not right...
#lambda = 0.3*lambda_max
#g = NormL21(lambda)

tol = 1e-4
Lf = 8

verb = 0
slv = PG(tol = tol, fast = true, gamma = 1/Lf, linesearch = false, verbose = verb)
slv = ZeroFPR(tol = tol, gamma = 1/Lf, linesearch = false, verbose = verb)
Y2 = zeros(Y)
lambdas = linspace(1e-5,0.08,100)*lambda_max
@time for lambda in lambdas
#	Y2 = zeros(Y)
	g = NormL1(lambda)
	Y2, slv = solve(Ladj, L, R_w, Conjugate(g), Y2, slv)
	show(slv)
end

Y2 = zeros(Y)
g = NormL1(lambdas[end])
Y2, slv = solve(Ladj, L, R_w, Conjugate(g), Y2, slv)
show(slv)

Y3 =-Ladj(Y2)+R_w

#ImageView.view(R,xy=["y","x"])
#ImageView.view(R_w,xy=["y","x"])
#ImageView.view(Y3,xy=["y","x"])

return
