using RegLS, AbstractOperators

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

A = [a; b]
#preprocess data
A .-= mean(A)
A ./= sqrt.(var(A,1))

L = 7  # inner layers
W1 = Variable(sqrt(2*L)*randn(2,L)) # first layer
S1 = Sigmoid((Na+Nb,L)) 
b1 = Variable(1)
W2 = Variable(sqrt(L*L)*randn(L,L)) # second layer
S2 = Sigmoid((Na+Nb,L)) 
b2 = Variable(1)
W3 = Variable(sqrt(L)*randn(L,1))   # third layer
S3 = Sigmoid((Na+Nb,1)) 
b3 = Variable(1)

nn = S3*((S2*((S1*(A*W1.+b1))*W2.+b2))*W3.+b3) #construct neural network operator

lambda = 3e-2
lambda1 = lambda/(2*L)
lambda2 = lambda/(L*L)
lambda3 = lambda/L

slv = ZeroFPR(tol = 1e-4) 
reg = lambda1*norm(W1,1)+lambda2*norm(W2,1)+lambda3*norm(W3,1)

slv = @minimize crossentropy(nn,bi)+reg with slv 
println(slv)

xx = linspace(-3,3,200)
yy = linspace(-3,3,200)
S1 = Sigmoid((L,)) 
S2 = Sigmoid((L,)) 
S3 = Sigmoid((1,)) 

regions = [(S3*((~W3)'*(S2*((~W2)'*(S1*((~W1)'*[xi;yi].+~b1)).+~b2)).+~b3))[1] 
	   for yi in yy, xi in xx]

using PyPlot
figure()
pcolormesh(xx,yy,regions, cmap = "Oranges" ) 
contour(xx,yy,regions, levels = 0.7*[maximum(regions)] ) 
plot(A[Na+1:end,1],A[Na+1:end,2], "ko")
plot(A[1:Na,1],A[1:Na,2], "r*")


