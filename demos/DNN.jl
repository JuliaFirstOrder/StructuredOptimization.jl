using RegLS

srand(123)
#construct dataset
Na,Nb = 700,700
rhoa = 5.
a = vcat([ rhoa.*[ cos(theta)  sin(theta)] for theta in linspace(0,2*pi,Na) ]...)   
a .+= 1.5.*randn(size(a))
b = [zeros(div(Nb,2))  zeros(div(Nb,2))]
rhob = 8.
b = [b; vcat([ rhob.*[ cos(theta)  sin(theta)] for theta in linspace(0,2*pi,div(Nb,2)) ]...)]   
b .+= randn(size(b))
bi = [ones(Bool,Na);zeros(Bool,Nb)]

A = [a' b']
#preprocess data
A .-= mean(A)
A ./= sqrt.(var(A,2))

L = 7  # inner layers
W1 = Variable(sqrt(2*L)*randn(L,2)) # first layer
S1 = Sigmoid((L,Na+Nb)) 
W2 = Variable(sqrt(L*L)*randn(L,L)) # second layer
S2 = Sigmoid((L,Na+Nb)) 
W3 = Variable(sqrt(L)*randn(1,L))   # third layer
S3 = Sigmoid((1,Na+Nb)) 

nn = S3*(W3*(S2*(W2*(S1*(W1*A))))) #construct neural network operator

lambda = 8e-5
lambda1 = lambda*sqrt(2*L)
lambda2 = lambda*sqrt(L*L)
lambda3 = lambda*sqrt(L)

slv = ZeroFPR(tol = 1e-4, maxit = 10000) 
reg = lambda1*norm(W1,1)+lambda2*norm(W2,1)+lambda3*norm(W3,1)

slv = @minimize crossentropy(nn,bi)+reg with slv 
println(slv)

xx = linspace(-3,3,200)
yy = linspace(-3,3,200)
S1 = Sigmoid((L,)) 
S2 = Sigmoid((L,)) 
S3 = Sigmoid((1,)) 

regions = [ (S3*(~W3*(S2*(~W2*(S1*(~W1*[xi;yi]))))))[1] for xi in xx, yi in yy]'

using PyPlot
figure()
pcolormesh(xx,yy,regions, cmap = "Oranges" ) 
contour(xx,yy,regions, levels = 0.7*[maximum(regions)] ) 
plot(A[1,Na+1:end],A[2,Na+1:end], "ko")
plot(A[1,1:Na],A[2,1:Na], "r*")


