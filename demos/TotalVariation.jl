using RegLS
using ProximalOperators
using Images
using ImageView
srand(123)

using TestImages
img = testimage("lena_gray")

R = convert(Array{Float64},img)
R_w = R+sqrt(0.006*norm(R[:],Inf))*randn(size(R))

slv = ZeroFPR
slv = FPG
tol = 1e-4
Lf  = 8
verb = 1
slv = slv(tol = tol, verbose = verb , gamma = 1/Lf, adaptive = false)
lambda = 0.07

X1 = Variable(zeros(R))
slv = @minimize ls(X1-R_w)+lambda*norm(variation(X1),1) with slv
println(slv)

X21 = Variable(zeros(R))
slv = @minimize ls(X21-R_w)+lambda*sum(norm(variation(X21)),2) with slv
println(slv)


imshow(R)
imshow(R_w)
imshow(~X1)
imshow(~X21)

