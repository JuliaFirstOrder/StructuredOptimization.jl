### Sparse signal reconstruction

srand(0) # for reproducibility

m, n, k = 1024, 4096, 160 # problem parameters
B = randn(m, n)
Q, R = qr(B')
A = Q' # measurement matrix
x_orig = sign(randn(n))
J = randperm(n)
x_orig[J[k+1:end]] = 0 # generate sparse signal
sigma = 1e-2
y = A*x_orig + sigma*randn(m) # add noise to measurement

using RegLS
lambda_max = norm(A'*y, Inf)
lambda = 0.01*lambda_max # regularization parameter
x_L1, info = solve(A, y, normL1(lambda), zeros(n))

x_L0c, info = solve(A, y, indBallL0(200), zeros(n))

return
