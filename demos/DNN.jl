using RegLS

srand(123)
Na,Nb = 700,700
rhoa = 5.
a = vcat([ rhoa.*[ cos(theta)  sin(theta)] for theta in linspace(0,2*pi,Na) ]...)   
a .+= 1.5.*randn(size(a))
b = [zeros(div(Nb,2))  zeros(div(Nb,2))]
rhob = 8.
b = [b; vcat([ rhob.*[ cos(theta)  sin(theta)] for theta in linspace(0,2*pi,div(Nb,2)) ]...)]   
b .+= randn(size(b))
bi = [ones(Bool,Na);zeros(Bool,Nb)]

A = [[a[:,1];b[:,1]] [a[:,2];b[:,2]]]
#preprocess data
A .-= mean(A)
A ./= sqrt.(var(A,1))

L = 7  #inner layers
W1 = Variable(sqrt(2*L)*randn(2,L))
S1 = Sigmoid((Na+Nb,L)) 
W2 = Variable(sqrt(L*L)*randn(L,L))
S2 = Sigmoid((Na+Nb,L)) 
W3 = Variable(sqrt(L)*randn(L,1))
S3 = Sigmoid((Na+Nb,1)) 

lambda = 8e-5
lambda1 = lambda*sqrt(2*L)
lambda2 = lambda*sqrt(L*L)
lambda3 = lambda*sqrt(L)

slv = ZeroFPR(tol = 1e-4, maxit = 20000)
reg = lambda1*norm(W1,1)+lambda2*norm(W2,1)+lambda3*norm(W3,1)

slv = @minimize crossentropy(  S3*(S2*((S1*(A*W1))*W2)*W3)  ,bi)+reg with slv 
println(slv)

xx = linspace(-3,3,200)
yy = linspace(-3,3,200)
S1 = Sigmoid((2,L)) 
S2 = Sigmoid((2,L)) 
S3 = Sigmoid((2,1)) 

using PyPlot
figure()
pcolormesh(xx,yy,[sum( S3*(S2*((S1*([xi yi]*~W1))*~W2)*~W3) )  for xi in xx, yi in yy]', cmap = "Oranges" ) 
contour(xx,yy,[sum( S3*(S2*((S1*([xi yi]*~W1))*~W2)*~W3) )  for xi in xx, yi in yy]')
plot(A[Na+1:end,1],A[Na+1:end,2], "ko")
plot(A[1:Na,1],A[1:Na,2], "r*")


