using RegLS

srand(123)
Na,Nb = 300,300
rho = 10.
a = vcat([ 10.*[ cos(theta)  sin(theta)] for theta in linspace(0,2*pi,Na) ]...)   
a .+= randn(size(a))
b = [zeros(div(Nb,2))  zeros(div(Nb,2))]
b = [b; vcat([ 20.*[ cos(theta)  sin(theta)] for theta in linspace(0,2*pi,div(Nb,2)) ]...)]   
b .+= randn(size(b))
bi = [ones(Bool,Na);zeros(Bool,Nb)]

gamma = 1.
L = 20
W1 = Variable(0.5*randn(2,L))
S1 = Sigmoid((Na+Nb,L),gamma) 
W2 = Variable(0.5*randn(L,L))
S2 = Sigmoid((Na+Nb,L),gamma) 
W3 = Variable(0.5*randn(L,1))
S3 = Sigmoid((Na+Nb,1)) 

A = [[a[:,1];b[:,1]] [a[:,2];b[:,2]]]
lambda = 1e-6
lambda1 = lambda*length(~W1)
lambda2 = lambda*length(~W2)
lambda3 = lambda*length(~W3)

slv = @minimize crossentropy(  S3*(S2*((S1*(A*W1))*W2)*W3)  ,bi)+lambda1*norm(W1,1)+lambda2*norm(W2,1)+lambda3*norm(W3,1) with ZeroFPR(tol = 1e-4, maxit = 2000) 

xx = linspace(-30,30,200)
yy = linspace(-30,30,200)
S1 = Sigmoid((2,L)) 
S2 = Sigmoid((2,L)) 
S3 = Sigmoid((2,1)) 

using PyPlot
figure()
pcolormesh(xx,yy,[sum( S3*(S2*((S1*([xi yi]*~W1))*~W2)*~W3) )  for xi in xx, yi in yy]', cmap = "Oranges" ) 
plot(a[:,1],a[:,2], "r*")
plot(b[:,1],b[:,2], "ko")


