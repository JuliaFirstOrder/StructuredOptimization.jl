x1,x2 = Variable(3),Variable(4)
X1,X2 = randn(3),randn(4)
X = Variable(3,4)
b1,b2 = randn(3),randn(4)
M1 = randn(4,3)
M2  = randn(3,4)

T = 0.1*ls(x1)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 0.1*ls(x1+M2*x2)
show(T)
y = RegLS.get_prox(T.f[1])
resx = affine(T)[1]([X1,X2])
println()

T = 0.1*ls(x1-M2*x2)+ls(x2)
show(T)
println()

T = norm(x1,0)
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])

println()
T = 0.1*norm(x1,1)
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = 0.1*norm(x1,2)
T = 0.1*norm(x1)
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = 0.1*norm(x1,Inf)
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1,0) <= 2
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = 3*norm(x1,1) <= 2
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1) <= 2
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = sum(norm(X),1)
RegLS.lambda(terms(T)[1])
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 0.5*sum(norm(X),2)
RegLS.lambda(terms(T)[1])
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1,Inf) <= 2
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1) == 2
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = x1 <= 2.0
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = 3.0 <= x1
show(T)
RegLS.lambda(terms(T)[1])
y = RegLS.get_prox(T.f[1])
println()

T = x1 >= randn(3)
RegLS.lambda(terms(T)[1])
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = x1 in [zeros(3), ones(3)]
RegLS.lambda(terms(T)[1])
show(T)
y = RegLS.get_prox(T.f[1])
println()

T =1.2* hingeloss(x1,randn(3))
RegLS.lambda(terms(T)[1])
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = rank(X) <= 2
RegLS.lambda(terms(T)[1])
show(T)
y = RegLS.get_prox(T.f[1])
println()


println("testing cost function constructors")

x,y,z,u = Variable(randn(4)), Variable(randn(4)), Variable(randn(2,2)), Variable(randn(4))
cf = ls(x)+ls(y)
@test length(variable(cf)) == 2
cf = ls(x+y)+ls(y)
@test length(variable(cf)) == 2
cf = ls(x+y)+ls(y)+ls(z)
@test length(variable(cf)) == 3
cf = ls(x)+ls(y)+ls(z)
@test length(variable(cf)) == 3
cf = 10*(4*ls(x)+ls(x+y)+2*ls(z))
@test length(variable(cf)) == 3
@test  cf.f[1].lambda == 40 && cf.f[2].lambda == 10 && cf.f[3].lambda == 20
show(cf)

println("\n testing gradient \n")
T = 0.1*ls(x1)
show(T)
resx, fx = RegLS.residual(T,X1)
fx = RegLS.residual!(resx, T, X1)
gradfi, fx = RegLS.gradient(T, resx)
@test norm(fx-(0.1/2*vecnorm(X1)^2))<1e-8
@test norm(gradfi[1]-0.1*resx[1])<1e-8
grad   = RegLS.At_mul_B(T, gradfi)
@test norm(grad-(0.1*(X1)))<1e-8
println()

T = 0.1*ls(x1)+ls(M1*x1)
show(T)
resx, fx = RegLS.residual(T,X1)
fx = RegLS.residual!(resx, T, X1)
gradfi, fx = RegLS.gradient(T, resx)
@test norm(fx-(0.1/2*vecnorm(X1)^2+1/2*vecnorm(M1*X1)^2))<1e-8
grad   = RegLS.At_mul_B(T, gradfi)
@test norm(grad-(0.1*(X1)+M1'*(M1*X1)))<1e-8
println()

T = 0.1*ls(x1+M2*x2-b1)+ls(M1*x1-x2-b2)
show(T)
resx, fx = RegLS.residual(T,[X1,X2])
fx = RegLS.residual!(resx, T, [X1,X2])
gradfi, fx = RegLS.gradient(T, resx)
@test norm(fx-(0.1/2*vecnorm(X1+M2*X2-b1)^2+1/2*vecnorm(M1*X1-X2-b2)^2))<1e-8
grad   = RegLS.At_mul_B(T, gradfi)
@test norm(grad-([0.1*(X1+M2*X2-b1)+M1'*(M1*X1-X2-b2), 0.1*M2'*(X1+M2*X2-b1)-(M1*X1-X2-b2)]))<1e-8


println("\n testing gradient Moreau \n")
gamma0  = 0.5
lambda = 0.1
T = lambda*norm(M1*x1-b2,1)
T = smooth(T, gamma0)
show(T)
resx, fx = RegLS.residual(T,X1)
fx2 = RegLS.cost(T, resx)
@test norm(resx[1]-(M1*X1-b2))<=1e-8
@test norm(fx2-fx)<1e-8
p, = prox(NormL1(lambda), resx[1],lambda*gamma0)
@test norm(fx-( lambda*norm(p,1)+1/(2*gamma0*lambda)*vecnorm(resx[1]-p)^2  ) )<1e-8
gradfi, fx3 = RegLS.gradient(T, resx)
@test norm(fx-fx3) <= 1e-8
grad   = RegLS.At_mul_B(T, gradfi)
grad2 = 1/(gamma0*lambda)*M1'*(resx[1]-p) 
@test norm(grad-grad2)<=1e-8

gamma0 = 0.5
T = lambda*norm(M1*x1+x2-b2,1)
T = smooth(T, gamma0)
show(T)
resx, fx = RegLS.residual(T,[X1,X2])
fx2 = RegLS.cost(T, resx)
@test norm(resx[1]-(M1*X1+X2-b2))<=1e-8
@test norm(fx2-fx)<1e-8
p, = prox(NormL1(lambda), resx[1],lambda*gamma0)
@test norm(fx-( lambda*norm(p,1)+1/(2*gamma0*lambda)*vecnorm(resx[1]-p)^2  ) )<1e-8
gradfi, fx3 = RegLS.gradient(T, resx)
@test norm(fx-fx3) <= 1e-8
grad   = RegLS.At_mul_B(T, gradfi)
gradX1 = 1/(gamma0*lambda)*M1'*(resx[1]-p)
gradX2 = 1/(gamma0*lambda)*(resx[1]-p)
@test norm(grad[1]-gradX1)<=1e-8
@test norm(grad[2]-gradX2)<=1e-8


println("\n test split \n")

x,y,z,u = Variable(randn(4)), Variable(randn(4)), Variable(randn(2,2)), Variable(randn(4))
cf = ls(x+y-randn(4))+ls(y+reshape(z,4)+u)+1*norm(randn(4).*u,1)+2*norm(z-randn(2,2),2)+3*norm(y,1)
@test length(variable(cf)) == 4
show(cf)

x_sorted = sort(variable(cf),by = object_id)

smooth, proximable, nonsmooth = RegLS.split(cf)
println("\n Smooth: \n")
show(smooth)
println("\n Proximable: \n")
show(proximable)
println("\n NonSmooth: \n")
show(nonsmooth)
println("\n  \n")

println("\n test merge prox \n")
p = RegLS.mergeProx(x_sorted, proximable)
show(x_sorted)
show(p.fs)

println("\n sort_and_expand \n")
smooth_exp = RegLS.sort_and_expand(x_sorted, smooth)
show(RegLS.affine(smooth))
show(RegLS.affine(smooth_exp))

out1 = RegLS.affine(smooth_exp)[1](~x_sorted)
out2 = RegLS.affine(smooth)[1]([~x,~y])
@test norm(out1-out2) <= 1e-8

out1 = RegLS.affine(smooth_exp)[2](~x_sorted)
out2 = RegLS.affine(smooth)[2](~[y,z,u])
@test norm(out1-out2) <= 1e-8
