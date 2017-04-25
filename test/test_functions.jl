println("\n testing functions \n")

x0 = randn(10)
lambda = 0.1
gamma = 1e-2
F = RegLS.LeastSquares(lambda)

grad, fx = RegLS.gradient(F,x0)

@test norm(grad-lambda*x0)<=1e-8
@test norm(fx-0.5*lambda*norm(x0)^2)<=1e-8

x1, fx1 = RegLS.gradstep(F,x0, gamma)
@test norm(x1-(x0 - gamma*lambda*x0))<=1e-8
@test norm(fx1- 0.5*lambda*norm(x0 - gamma*lambda*x0)^2 )<=1e-8

@test norm(F(x1)-fx1) < 1e-8

###########-------------###########

x0 = randn(10)
lambda1 = maximum(abs(x0))*0.3
lambda  = pi 
gamma  = 1e-2
F = RegLS.MoreauEnvelope(lambda,NormL1(lambda1))

grad, fx = RegLS.gradient(F,x0)

pr, = prox(NormL1(lambda1), x0, lambda)

@test norm(grad - (  1/lambda*(x0-pr)   )   ) <=1e-8
@test norm(fx - ( lambda1*norm(pr,1)+1/(2*lambda)*norm(x0-pr)^2   )   ) <=1e-8

x1, fx1 = RegLS.gradstep(F,x0, gamma)
@test norm(x1 - ( x0 - gamma/lambda*(x0-pr)   )   ) <=1e-8
pr1, = prox(NormL1(lambda1), x1, lambda)
@test norm(fx1 - ( lambda1*norm(pr1,1)+1/(2*lambda)*norm(x1-pr1)^2   )   ) <=1e-8

@test norm(F(x1)-fx1) < 1e-8
