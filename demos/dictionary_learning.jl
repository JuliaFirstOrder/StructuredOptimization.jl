
using RegLS
using ProximalOperators
#using Images, ImageView
srand(123)

#using TestImages
#img = testimage("lena_gray")


n,m = 10,20
l = 25
Y = randn(n,m)
D = randn(n,l)
X = randn(l,m)

DX = [D,X]
L = DX -> DX[1]*DX[2]
Ladj = Y -> [(DX[2]*Y')',DX[1]'*Y]
norm(vecdot(L(DX),Y)-vecdot(DX[1],Ladj(Y)[1])) #verify adjoint operator
norm(vecdot(L(DX),Y)-vecdot(DX[2],Ladj(Y)[2])) #verify adjoint operator


D = randn(n,l)
for i = 1:l
	D[:,i] = D[:,i]./norm(D[:,i])
end
X = full(sprandn(l,m,0.1))
Y = D*X+0.001*randn(n,m) 

gg = SlicedSeparableSum(repmat([IndSphereL2(1.)],l),[[i] for i = 1:l ],2)
g = SeparableSum([gg,IndBallL0(countnz(X)+10) ])
DX = [ones(D),ones(X)]
DX, =  solve!(L, Ladj, Y, g, DX, ZeroFPR(verbose = 1, tol = 1e-7, gamma = 1e-7, linesearch = false))
