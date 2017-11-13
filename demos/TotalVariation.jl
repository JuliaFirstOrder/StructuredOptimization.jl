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
#slv = FPG
tol = 1e-4
Lf  = 8
verb = 1
slv = slv(tol = tol, verbose = verb, gamma = 1/Lf, adaptive = false)
lambda = 0.07

L = Variation(Float64,size(R))

Y1 = Variable(size(L,1)...)
slv = @minimize ls(L'*Y1-R_w)+lambda*conj(norm(Y1,1)) with slv
println(slv)
R1 = -(L'*(~Y1)-R_w)

Y1 = Variable(size(L,1)...)
slv = @minimize ls(L'*Y1-R_w)+lambda*conj(norm(Y1,2,1,2)) with slv
println(slv)
R21 = -(L'*(~Y1)-R_w)

imshow(R)
imshow(R_w)
imshow(R1)
imshow(R21)

