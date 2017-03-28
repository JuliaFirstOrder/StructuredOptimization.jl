x = OptVar(3)
X = randn(3)

T = 0.1*ls(x)
show(T)
y = RegLS.get_prox(T.f[1])
@test norm(T(X)-0.1/2*vecnorm(X)^2)<1e-8
println()

x,y = OptVar(3),OptVar(3)
X,Y = randn(3),randn(3)
M = randn(3,3)

T = 0.1*ls(M*x+y)
show(T)
y = RegLS.get_prox(T.f[1])
@test norm(T([X,Y])-0.1/2*vecnorm(M*X+Y)^2)<1e-8
resx = affine(T)[1]([X,Y])
@test norm(T([X,Y],false)-0.1/2*vecnorm(M*X+Y)^2)<1e-8
println()

#T = norm(x,0)
#show(T)
#y = RegLS.get_prox(T.f[1])
#
#println()
#T = 0.1*norm(x,1)
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = 0.1*norm(x,2)
#T = 0.1*norm(x)
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = 0.1*norm(x,Inf)
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = norm(x,0) <= 2
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = 3*norm(x,1) <= 2
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = norm(x) <= 2
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = norm(x,Inf) <= 2
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = norm(x) == 2
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = x <= 2.0
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = 3.0 <= x
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = x >= randn(3)
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T = x in [zeros(3), ones(3)]
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
#T =1.2* hingeloss(x,randn(3))
#show(T)
#y = RegLS.get_prox(T.f[1])
#println()
#
