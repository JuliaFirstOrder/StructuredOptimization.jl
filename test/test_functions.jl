x1,x2 = OptVar(3),OptVar(4)
X1,X2 = randn(3),randn(4)
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

T = 0.1*ls(x1)
show(T)
resx, fx = RegLS.evaluate(T,X1)
@test norm(fx-(0.1/2*vecnorm(X1)^2))<1e-8
fx  = RegLS.evaluate!(resx, T, X1)
fx2 = RegLS.cost(T, resx)
@test norm(fx2-fx)<1e-8
@test norm(fx-(0.1/2*vecnorm(X1)^2))<1e-8
grad = RegLS.gradient(T, resx)
@test norm(grad-(0.1*(X1)))<1e-8
println()

T = 0.1*ls(x1)+ls(M1*x1)
show(T)
resx, fx = RegLS.evaluate(T,X1)
@test norm(fx-(0.1/2*vecnorm(X1)^2+1/2*vecnorm(M1*X1)^2))<1e-8
fx = RegLS.evaluate!(resx, T, X1)
fx2 = RegLS.cost(T, resx)
@test norm(fx2-fx)<1e-8
@test norm(fx-(0.1/2*vecnorm(X1)^2+1/2*vecnorm(M1*X1)^2))<1e-8
grad = RegLS.gradient(T, resx)
@test norm(grad-(0.1*(X1)+M1'*(M1*X1)))<1e-8
println()

T = 0.1*ls(x1+M2*x2-b1)+ls(M1*x1-x2-b2)
show(T)
resx, fx = RegLS.evaluate(T,[X1,X2])
@test norm(fx-(0.1/2*vecnorm(X1+M2*X2-b1)^2+1/2*vecnorm(M1*X1-X2-b2)^2))<1e-8
fx = RegLS.evaluate!(resx, T, [X1,X2])
fx2 = RegLS.cost(T, resx)
@test norm(fx2-fx)<1e-8
@test norm(fx-(0.1/2*vecnorm(X1+M2*X2-b1)^2+1/2*vecnorm(M1*X1-X2-b2)^2))<1e-8
grad = RegLS.gradient(T, resx)
@test norm(grad-([0.1*(X1+M2*X2-b1)+M1'*(M1*X1-X2-b2), 0.1*M2'*(X1+M2*X2-b1)-(M1*X1-X2-b2)]))<1e-8
println()

T = norm(x1,0)
show(T)
y = RegLS.get_prox(T.f[1])

println()
T = 0.1*norm(x1,1)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 0.1*norm(x1,2)
T = 0.1*norm(x1)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 0.1*norm(x1,Inf)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1,0) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 3*norm(x1,1) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1,Inf) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x1) == 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = x1 <= 2.0
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 3.0 <= x1
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = x1 >= randn(3)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = x1 in [zeros(3), ones(3)]
show(T)
y = RegLS.get_prox(T.f[1])
println()

T =1.2* hingeloss(x1,randn(3))
show(T)
y = RegLS.get_prox(T.f[1])
println()

