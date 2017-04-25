println("\n testing functions \n")

x0 = randn(10)
lambda = 0.1
gamma = 1e-2
f = SqrNormL2(lambda)

grad, fx = RegLS.gradient(f, x0)

@test norm(grad - lambda*x0) <= 1e-8
@test norm(fx - 0.5*lambda*norm(x0)^2) <= 1e-8

x1, fx0 = RegLS.gradstep(f, x0, gamma)

@test norm(x1 - (x0 - gamma*lambda*x0)) <= 1e-8
@test norm(fx0 - 0.5*lambda*norm(x0)^2 ) <= 1e-8
@test abs(f(x0) - fx0) <= 1e-8

###########-------------###########

x0 = randn(10)
mu = maximum(abs(x0))*0.3
lambda = pi
gamma = 1e-2
f = RegLS.MoreauEnvelope(lambda, NormL1(mu))

grad, fx = RegLS.gradient(f, x0)

pr, = prox(NormL1(mu), x0, lambda)

@test norm(grad - (  1.0/lambda*(x0-pr)   )   ) <= 1e-8
@test norm(fx - ( mu*norm(pr,1)+1.0/(2.0*lambda)*norm(x0-pr)^2   )   ) <= 1e-8

x1, fx0 = RegLS.gradstep(f, x0, gamma)

@test norm(x1 - ( x0 - gamma/lambda*(x0-pr)   )   ) <= 1e-8
@test norm(fx0 - ( mu*norm(pr,1)+1.0/(2.0*lambda)*norm(x0-pr)^2   )   ) <= 1e-8
@test abs(f(x0) - fx0) <= 1e-8
