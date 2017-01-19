
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
R_w = R+sqrt(0.001*norm(R[:],Inf))*randn(Nx,Ny)

M = 2^3
D = zeros(M,M,M^2)
counter = 1
for p = 0:M-1,q = 0:M-1
	alphap = p==0 ? 1/sqrt(M) : sqrt(2/M)
	alphaq = q==0 ? 1/sqrt(M) : sqrt(2/M)
	D[:,:,counter] = (
	[alphap*alphaq*cos(pi*(2*m+1)*p/(2*M))*cos(pi*(2*n+1)*q/(2*M)) for m=0:M-1,n=0:M-1])
	counter +=1
end
D = reshape(D,M^2,M^2)


D2 = [D randn(M^2,3*M^2)]
#D2 = [D ones(M^2,3*M^2)]
#D2 = repmat(D,1,4)
#D2 =  ones(M^2,4*M^2)
D2 = copy(D) 
Na = size(D2,2)
Np = round(Int64,Nx*Ny/M^2) #number of patches
Y = zeros(M^2,Np)

c = 0
c2 = 0
for i = 1:Np
	Y[:,i] = R_w[1+c:c+M,1+c2:c2+M][:] 
	c += M 
	if c == 256
		c = 0
		c2+=M
	end
end

X = zeros(Na,Np)
DX = [D2,X]
L = DX -> DX[1]*DX[2]
Ladj = Y -> [(DX[2]*Y')',DX[1]'*Y]
#norm(vecdot(L(DX),Y)-vecdot(DX[1],Ladj(Y)[1])) #verify adjoint operator
#norm(vecdot(L(DX),Y)-vecdot(DX[2],Ladj(Y)[2])) #verify adjoint operator

#gg  = SlicedSeparableSum(repmat([IndSphereL2(1.)],Na),[[i] for i = 1:Na ],2)
#ggg = SlicedSeparableSum(repmat([IndBallL0(50)],Np),[[i] for i = 1:Np ],2)
##ggg = NormL0(1e-4)
#g = SeparableSum([gg,ggg])
#
#DX, = solve!(L, Ladj, Y, g, DX, ZeroFPR(verbose = 2, tol = 1e-8, gamma = 1e-6,  linesearch = false, maxit = 1000))
#
#Y2 = L(DX)
#
#D = reshape(D,M,M,M^2)
#D2 = reshape(DX[1],M,M,Na)
#
#D_view = zeros(M^2,M^2)
#D2_view = zeros(M^2,Na)
#c = 0
#c2 = 0
#for i = 1:M^2
#	D_view[1+c:c+M,1+c2:c2+M] = D[:,:,i]
#	c += M 
#	if c == M^2
#		c = 0
#		c2+=M
#	end
#end
#
#c = 0
#c2 = 0
#for i = 1:Na
#	D2_view[1+c:c+M,1+c2:c2+M] = D2[:,:,i]
#	c += M 
#	if c == M^2
#		c = 0
#		c2+=M
#	end
#end




#lasso on DCT dictionary
D = reshape(D,M^2,M^2)
#D2 = reshape(D2,M^2,Na)
#D = D2

X = zeros(size(D,2),Np)
L = X -> D*X
Ladj = Y -> D'*Y

lambda_max = vecnorm(Ladj(Y), Inf)
#g = SlicedSeparableSum(repmat([IndBallL0(15)],Np),[[i] for i = 1:Np ],2)
g = NormL1(lambda_max*0.005)
X, = solve(L, Ladj, Y, g, X, ZeroFPR(verbose = 2))

Y2 = D*X

R2 = zeros(R)
c = 0
c2 = 0
for i = 1:Np
	R2[1+c:c+M,1+c2:c2+M] = reshape(Y2[:,i],M,M)
	c += M 
	if c == 256
		c = 0
		c2+=M
	end
end


ImageView.view(R,xy=["y","x"])
ImageView.view(R_w,xy=["y","x"])
ImageView.view(R2,xy=["y","x"])








