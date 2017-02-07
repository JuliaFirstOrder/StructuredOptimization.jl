using RegLS
using ProximalOperators
using Images
using ImageView
srand(123)

using TestImages
img = testimage("lena_gray")

R = convert(Array{Float64},img)

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

tol = 1e-4
Lf = 8

verb = 1
slv1 = PG(tol = tol, fast = true, gamma = 1/Lf, linesearch = false, verbose = verb)
slv2 = ZeroFPR(tol = tol, gamma = 1/Lf, linesearch = false, verbose = verb)

Y2 = zeros(Y)

g = Regularize(NormL1(lambda),1-lambda)
#g = NormL1(lambda)
@time Y3,slv1 = solve(R_w, g, L, Ladj, Y2, slv1)
println("prox eval:$(slv1.cnt_prox), matvec: $(slv1.cnt_matvec)")

@time Y3,slv2 = solve(R_w, g, L, Ladj, Y2, slv2)
println("prox eval:$(slv2.cnt_prox), matvec: $(slv2.cnt_matvec)")

imshow(R)
imshow(R_w)
imshow(Y3)

return
