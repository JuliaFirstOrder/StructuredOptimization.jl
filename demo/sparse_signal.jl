using RegLS

# for reproducibility
srand(0)

# problem parameters
n = 4096
m = 1024
k = 160

# measurement matrix
B = randn(m, n)
Q, R = qr(B')
A = Q'

# generate signal and noise
x_orig = sign(randn(n))
J = randperm(n)
x_orig[J[k+1:end]] = 0
sigma = 1e-2
y = A*x_orig + sigma*randn(m)

# solve problem
lambda_max = norm(A'*y, Inf)
lambda = 0.01*lambda_max
x_L1, it = solve(A, y, normL1(lambda), zeros(n))
x_L0c, it = solve(A, y, indBallL0(k), zeros(n))

return
