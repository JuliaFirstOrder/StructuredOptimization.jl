println("\nTesting cost terms\n")

# Simple Terms

x = Variable(randn(10))
X = Variable(randn(3,4))
A = randn(4,10)
b = randn(4)

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
@test cf.f(~X) == sum(  sqrt.(sum((~X).^2, dims=1 )) ) 

cf = 4*norm(X,2,1,2)
@test cf.lambda - 4 == 0
@test cf.f(~X) == sum(  sqrt.(sum((~X).^2, dims=2 )) ) 

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

cf = 2*rank(X) <= 6
@test cf.lambda == 1
@test cf.f(~X) == (IndBallRank(3))(~X)

cf = rank(X)
@test_throws MethodError cf.f(~X)

cf = norm(X,*)
U, S, V = svd(~X)
@test cf.lambda == 1
@test cf.f(~X) == sum(S)

cf = rank(X)
@test_throws MethodError cf.f(~X)

y = randn(size(~x))
cf = hingeloss(x,y)
@test cf.lambda == 1
@test cf.f(~x) == (HingeLoss(y))(~x)

y = randn(size(~x))
cf = sqrhingeloss(x,y)
@test cf.lambda == 1
@test cf.f(~x) == (SqrHingeLoss(y))(~x)

y = randn(size(~x))
cf = logisticloss(x,y)
@test cf.lambda == 1
@test cf.f(~x) == (LogisticLoss(y))(~x)

xp = Variable(rand(10)) 
bp = rand(Float64, size(~xp))
cf = crossentropy(xp,bp)
@test cf.lambda == 1
@test cf.f(~xp) == (CrossEntropy(bp))(~xp)

cf = logbarrier(x)
@test cf.lambda == 1
@test cf.f(~x) == (LogBarrier(1.0))(~x)

cf = maximum(x)
@test cf.lambda == 1
@test cf.f(~x) == (Maximum(1.0))(~x)

cf = sumpositive(x)
@test cf.lambda == 1
@test cf.f(~x) == (SumPositive())(~x)

a = 1.
cf = huberloss(x,a)
@test cf.lambda == 1
@test cf.f(~x) == (HuberLoss(a))(~x)

a = randn(size(x))
cf = dot(a,x)
@test cf.lambda == 1
@test cf.f(~x) == (Linear(a))(~x)

#IndBinary
lu = (-1.0,randn(length(~x)))
cf = x == lu
@test cf.lambda == 1
@test cf.f(~x) == (IndBinary(lu...))(~x)

# IndAffine (not working in julia < 1.1)
if VERSION.major >= 1 && VERSION.minor >= 1
  cf = A*x-b == 0
  @test cf.lambda == 1
  @test cf.f(~x) == (IndAffine(A,b))(~x)

  cf = (A*x == b)
  @test cf.lambda == 1
  @test cf.f(~x) == (IndAffine(A,-b))(~x)
end

cf = 2*norm(x,1)
ccf = conj(cf)
@test ccf.A == cf.A
@test ccf.f == Conjugate(Postcompose(NormL1(),2))
@test_throws ErrorException conj(norm(randn(2,10)*x,1))

cf = 2*norm(x,1)
ccf = smooth(cf,2.0)
@test ccf.A == cf.A
@test ccf.f(~x) == MoreauEnvelope(Postcompose(NormL1(),2),2.0)(~x)

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
@test norm(affine(cf[1])*(~x) - (A*(~x)-b)) < 1e-12
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
@test StructuredOptimization.is_smooth(cf) == false
@test StructuredOptimization.is_smooth(cf^2) == true

cf = norm(w + z)^2
@test StructuredOptimization.is_smooth(cf) == true
@test StructuredOptimization.is_AcA_diagonal(cf) == false

cf = norm(x, 1) + norm(y, 2)
@test StructuredOptimization.is_smooth.(cf) == (false,false)
@test StructuredOptimization.is_AcA_diagonal.(cf) == (true,true)
