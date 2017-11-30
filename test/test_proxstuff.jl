
##TOMOVE to ProximalGradient
#
#x = randn(5)
#b = randn(5)
#mu = rand()
#
#f = SqrHingeLoss(b,mu)
#grad, f2 = gradient(f,x)
#
#gradfd = zeros(grad)
#for i in eachindex(gradfd)
#	delta = zeros(grad)
#	delta[i] = sqrt(eps())
#	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
#end
#
#@test f(x) == f2
#@test norm(grad-gradfd)/norm(grad) < 1e-6
#
#x = rand(5)
#b = [zeros(3);ones(2)]
#f = CrossEntropy(b)
#grad, f2 = gradient(f,x)
#
#gradfd = zeros(grad)
#for i in eachindex(gradfd)
#	delta = zeros(grad)
#	delta[i] = sqrt(eps())
#	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
#end
#
#@test norm(grad-gradfd)/norm(grad) < 1e-6
#@test norm(f(x) - f2) < 1e-7
#
#x = rand(1,5)
#b = [zeros(3);ones(2)]
#f = CrossEntropy(b)
#grad, f2 = gradient(f,x)
#
#gradfd = zeros(grad)
#for i in eachindex(gradfd)
#	delta = zeros(grad)
#	delta[i] = sqrt(eps())
#	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
#end
#
#@test norm(grad-gradfd)/norm(grad) < 1e-6
#@test norm(f(x) - f2) < 1e-7
#
#x = rand(5)
#b = [zeros(Bool,3);ones(Bool,2)]
#f = CrossEntropy(b)
#grad, f2 = gradient(f,x)
#
#gradfd = zeros(grad)
#for i in eachindex(gradfd)
#	delta = zeros(grad)
#	delta[i] = sqrt(eps())
#	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
#end
#
#@test norm(f(x) - f2) < 1e-7
#@test norm(grad-gradfd)/norm(grad) < 1e-6
#
#x = rand(1,5)
#b = [zeros(Bool,3);ones(Bool,2)]
#f = CrossEntropy(b)
#grad, f2 = gradient(f,x)
#
#gradfd = zeros(grad)
#for i in eachindex(gradfd)
#	delta = zeros(grad)
#	delta[i] = sqrt(eps())
#	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
#end
#
#@test norm(f(x) - f2) < 1e-7
#@test norm(grad-gradfd)/norm(grad) < 1e-6


# Testing precomposition by nonlinear operator

b = randn(10)
g = ProximalOperators.Translate(ProximalOperators.SqrNormL2(3.0), -b)
G = AbstractOperators.Sigmoid((10,), 1.0)
f = RegLS.PrecomposeNonlinear(g, G)

x = randn(10)

grad_f_x, f_x = gradient(f, x)

@test size(grad_f_x) == size(x)
@test abs(f_x - 3.0/2 * vecnorm(1.0 ./ (1.0 + exp.(-x)) - b)^2) <= 1e-10
expx = exp.(x)
expmx = 1.0./expx
grad_f_x_ref = 3.0 * ( expx ./ (1 + expx).^2 ) .* (1.0 ./ (1.0 .+ expmx) - b)
@test vecnorm(grad_f_x - grad_f_x_ref) <= 1e-10
