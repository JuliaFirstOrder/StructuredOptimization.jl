# Testing precomposition by nonlinear operator

b = randn(10)
g = ProximalOperators.Translate(ProximalOperators.SqrNormL2(3.0), -b)
G = AbstractOperators.Sigmoid((10,), 1.0)
f = StructuredOptimization.PrecomposeNonlinear(g, G)

x = randn(10)

grad_f_x, f_x = gradient(f, x)

@test size(grad_f_x) == size(x)
@test abs(f_x - 3.0/2 * vecnorm(1.0 ./ (1.0 + exp.(-x)) - b)^2) <= 1e-10
expx = exp.(x)
expmx = 1.0./expx
grad_f_x_ref = 3.0 * ( expx ./ (1 + expx).^2 ) .* (1.0 ./ (1.0 .+ expmx) - b)
@test vecnorm(grad_f_x - grad_f_x_ref) <= 1e-10
