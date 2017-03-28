x = OptVar(3)
X = randn(3)

T = 0.1*ls(x)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x,0)
show(T)
y = RegLS.get_prox(T.f[1])

println()
T = 0.1*norm(x,1)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 0.1*norm(x,2)
T = 0.1*norm(x)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 0.1*norm(x,Inf)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x,0) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 3*norm(x,1) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x,Inf) <= 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = norm(x) == 2
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = x <= 2.0
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = 3.0 <= x
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = x >= randn(3)
show(T)
y = RegLS.get_prox(T.f[1])
println()

T = x in [zeros(3), ones(3)]
show(T)
y = RegLS.get_prox(T.f[1])
println()

T =1.2* hingeloss(x,randn(3))
show(T)
y = RegLS.get_prox(T.f[1])
println()

