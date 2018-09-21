# Testing precomposition by nonlinear operator

b = randn(10)
g = ProximalOperators.SqrNormL2(3.0)
G = AffineAdd(AbstractOperators.Sigmoid((10,), 1.0),b,false)
f = StructuredOptimization.PrecomposeNonlinear(g, G)

x = randn(10)

grad_f_x, f_x = ProximalOperators.gradient(f, x)

@test size(grad_f_x) == size(x)
@test abs(f_x - 3.0/2 * norm(1.0 ./ (1.0 .+ exp.(-x)) - b)^2) <= 1e-10
expx = exp.(x)
expmx = 1.0./expx
grad_f_x_ref = 3.0 * ( expx ./ (1 .+ expx).^2 ) .* (1.0 ./ (1.0 .+ expmx) - b)
@test norm(grad_f_x - grad_f_x_ref) <= 1e-10

## with compose
#with vectors
l,m1,m2,n1,n2 = 2,3,4,5,6
x = (randn(m1,m2),randn(n1,n2))
A = randn(l,m1)
B = randn(m2,n1)
r = randn(l,n2)

b = randn(l,n2)
G = AffineAdd(NonLinearCompose( MatrixOp(A,m2), MatrixOp(B,n2) ),b,false)

g = ProximalOperators.SqrNormL2(3.0)
f = StructuredOptimization.PrecomposeNonlinear(g, G)

x = (randn(m1,m2),randn(n1,n2))

grad_f_x, f_x = ProximalOperators.gradient(f, x)

r = G*x 
grad_f_x2, f_x2 = ProximalOperators.gradient(g, r)
grad_f_x2 = jacobian(G,x)'*grad_f_x2

@test norm(f_x-f_x2) < 1e-8
@test blockvecnorm(grad_f_x2.-grad_f_x2) < 1e-8
