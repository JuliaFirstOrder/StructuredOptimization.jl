@printf("\nTesting cost terms\n")

# Simple Terms

x = Variable(randn(10))
X = Variable(randn(3,4))
A = randn(4,10)

cf = norm(x, 0)
@test cf.lambda == 1
@test cf.f(~x) == norm(~x,0)

cf = 3*norm(x, 0)
@test cf.lambda == 3
@test cf.f(~x) == norm(~x,0)

cf = norm(x, 0) <= 3
@test cf.lambda == 1
@test cf.f(~x) == (IndBallL0(3))(~x)

cf = norm(x, 1)
@test cf.lambda == 1
@test cf.f(~x) == (NormL1())(~x)

cf = norm(x, 1) <= 1.5
@test cf.lambda == 1
@test cf.f(~x) == (IndBallL1(1.5))(~x)

cf = 10*norm(x, 1) <= 1.5
@test cf.lambda == 1
@test cf.f(~x) == (IndBallL1(1.5/10))(~x)

cf = norm(x)
@test cf.lambda == 1
@test cf.f(~x) == norm(~x)

cf = pi*norm(x,2)
@test cf.lambda - pi == 0
@test cf.f(~x) == norm(~x)

cf = 3*norm(X,2,1)
@test cf.lambda - 3 == 0
@test cf.f(~X) == sum(sqrt.(sum((~X).^2,1))) 

cf = 4*norm(X,2,1,2)
@test cf.lambda - 4 == 0
@test cf.f(~X) == sum(sqrt.(sum((~X).^2,2))) 

@test_throws ErrorException 4*norm(X,1,2)

cf = norm(x, 2) <= 2.3
@test cf.lambda == 1
@test cf.f(~x) == (IndBallL2(2.3))(~x)

cf = norm(x, 2) == 2.3
@test cf.lambda == 1
@test cf.f(~x) == (IndSphereL2(2.3))(~x)

cf = norm(x, Inf)
@test cf.lambda == 1
@test cf.f(~x) == norm(~x,Inf)

cf = norm(x, Inf) <= 5.0
@test cf.lambda == 1
@test cf.f(~x) == (IndBallLinf(5.0))(~x)

cf = x <= 3.0
@test cf.lambda == 1
@test cf.f(~x) == (IndBox(-Inf, 3.0))(~x)

cf = 3.0 <= x
@test cf.lambda == 1
@test cf.f(~x) == (IndBox(3.0, Inf))(~x)

cf = x >= 1.0
@test cf.lambda == 1
@test cf.f(~x) == (IndBox(1.0, Inf))(~x)

cf = 1.0 >= x
@test cf.lambda == 1
@test cf.f(~x) == (IndBox(-Inf,1.0))(~x)

cf = x in [-5.0, 5.0]
@test cf.lambda == 1
@test cf.f(~x) == (IndBox(-5.0, 5.0))(~x)

cf = norm(x, 2)^2
@test cf.lambda == 1
@test cf.f(~x) == norm(~x)^2

cf = 0.5*norm(x, 2)^2
@test cf.lambda == 0.5
@test cf.f(~x) == norm(~x)^2

cf = 7*(0.5*norm(x, 2))^2
@test cf.lambda == 7*0.25
@test cf.f(~x) == norm(~x)^2

X = Variable(10, 10)
cf = 2*rank(X) <= 6
@test cf.lambda == 1
@test cf.f(~X) == (IndBallRank(3))(~X)

cf = rank(X)
@test_throws MethodError cf.f(~X)

b = randn(size(~x))
cf = hingeloss(x,b)
@test cf.lambda == 1
@test cf.f(~x) == (HingeLoss(b))(~x)

b = randn(size(~x))
cf = sqrhingeloss(x,b)
@test cf.lambda == 1
@test cf.f(~x) == (SqrHingeLoss(b))(~x)

x = Variable(rand(10))
b = rand(size(~x))
cf = crossentropy(x,b)
@test cf.lambda == 1
@test cf.f(~x) == (CrossEntropy(b))(~x)

x = Variable(randn(10))
a = 1.
cf = logbarrier(x,a)
@test cf.lambda == 1
@test cf.f(~x) == (LogBarrier(a))(~x)

a = 1.
cf = huberloss(x,a)
@test cf.lambda == 1
@test cf.f(~x) == (HuberLoss(a))(~x)

cf = 2*norm(x,1)
ccf = conj(cf)
@test ccf.A == cf.A
@test ccf.f == Conjugate(Postcompose(NormL1(),2))
@test_throws ErrorException conj(norm(randn(2,10)*x,1))

# Summing terms

x = Variable(10)
cf = ls(x) + 10*norm(x, 1)
@test cf[1].lambda == 1
@test cf[1].f(~x) == 0.5*norm(~x)^2
@test cf[2].lambda == 10
@test cf[2].f(~x) == norm(~x,1)

x = Variable(10)
cf = () #empty cost function 
cf += 10*norm(x, 1)
@test length(cf) == 1
@test cf[1].lambda == 10
@test cf[1].f(~x) == 10*norm(~x,1)

x = Variable(10)
cf = () #empty cost function 
cf += ls(x) + 10*norm(x, 1)
@test cf[1].lambda == 1
@test cf[1].f(~x) == 0.5*norm(~x)^2
@test cf[2].lambda == 10
@test cf[2].f(~x) == norm(~x,1)

# More complex situations

x = Variable(10)
A = randn(5, 10)
y = Variable(7)
B = randn(5, 7)
b = randn(5)

cf = ls(A*x - b) + norm(x, 1)
@test cf[1].lambda == 1
@test cf[1].f(~x) == 0.5*norm(~x)^2
@test displacement(cf[1]) == -b
@test operator(cf[1])*(~x) == A*(~x)
@test cf[2].lambda == 1
@test cf[2].f(~x) == norm(~x,1)

cf = ls(A*x - B*y + b) + norm(y, 1) + 5*norm(y, 2)
@test cf[1].lambda == 1
@test cf[1].f(~x) == 0.5*norm(~x)^2
@test cf[2].lambda == 1
@test cf[2].f(~x) == norm(~x,1)
@test cf[3].lambda == 5
@test cf[3].f(~x) == norm(~x,2)

cf = 10*(ls(A*x - B*y + b) + norm(y, 1) + 5*norm(y, 2))
@test cf[1].lambda == 10
@test cf[1].f(~x) == 0.5*norm(~x)^2
@test cf[2].lambda == 10
@test cf[2].f(~x) == norm(~x,1)
@test cf[3].lambda == 50
@test cf[3].f(~x) == norm(~x,2)

cf = 0.5*norm(A*x - B*y + b, 2)^2 + norm(x, 1) + norm(y, 2)
@test cf[1].lambda == 0.5
@test cf[1].f(~x) == norm(~x)^2
@test cf[2].lambda == 1
@test cf[2].f(~x) == norm(~x,1)
@test cf[3].lambda == 1
@test cf[3].f(~x) == norm(~x,2)

# Properties
A = randn(5, 10)
u = Variable(5)
w = Variable(5)
z = Variable(5)

cf = norm(A*x + z)
@test RegLS.is_smooth(cf) == false
@test RegLS.is_smooth(cf^2) == true

cf = norm(w + z)^2
@test RegLS.is_smooth(cf) == true
@test RegLS.is_AcA_diagonal(cf) == false

cf = norm(x, 1) + norm(y, 2)
@test RegLS.is_smooth.(cf) == (false,false)
@test RegLS.is_AcA_diagonal.(cf) == (true,true)
