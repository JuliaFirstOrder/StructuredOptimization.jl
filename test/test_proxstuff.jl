
#TOMOVE to ProximalGradient

x = randn(5)
b = randn(5)
mu = rand()

f = SqrHingeLoss(b,mu)
grad, f2 = gradient(f,x)

gradfd = zeros(grad)
for i in eachindex(gradfd)
	delta = zeros(grad)
	delta[i] = sqrt(eps())
	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
end

@test f(x) == f2
@test norm(grad-gradfd)/norm(grad) < 1e-6

x = rand(5)
b = [zeros(3);ones(2)]
f = CrossEntropy(b)
grad, f2 = gradient(f,x)

gradfd = zeros(grad)
for i in eachindex(gradfd)
	delta = zeros(grad)
	delta[i] = sqrt(eps())
	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
end

@test f(x) == f2
@test norm(grad-gradfd)/norm(grad) < 1e-6

x = rand(5)
b = [zeros(Bool,3);ones(Bool,2)]
f = CrossEntropy(b)
grad, f2 = gradient(f,x)

gradfd = zeros(grad)
for i in eachindex(gradfd)
	delta = zeros(grad)
	delta[i] = sqrt(eps())
	gradfd[i] = (f(x+delta)-f(x))/sqrt(eps())
end

@test f(x) == f2
@test norm(grad-gradfd)/norm(grad) < 1e-6
