n = 2100 # number of features
m = 130 # number of data points

w = randn(n)

A = randn(m, n)
btrue = sign(A*w)

b = sign(btrue + sqrt(0.1)*randn(m))
mu = 1.0

# solve dual problem (it's a QP)
BA = b.*A;
Q = BA*BA';
q = ones(size(A, 1));
sol = quadprog(-q, Q, eye(m), '<', Inf, 0.0, mu, IpoptSolver(print_level=0))
x_qp = BA'*sol.sol;

println("Solving random SVM problem: default solver/options")
x, slv = solve(zeros(n), HingeLoss(b, mu), A, zeros(m))
@test norm(x-x_qp, Inf)/norm(x_qp, Inf) <= 1e-5

println("Solving random SVM problem: random initial point")
x, slv = solve(zeros(n), HingeLoss(b, mu), A, randn(m))
@test norm(x-x_qp, Inf)/norm(x_qp, Inf) <= 1e-5

println("Solving random SVM problem: proximal gradient (quiet)")
x, slv = solve(zeros(n), HingeLoss(b, mu), A, zeros(m), PG(verbose = 0))
@test norm(x-x_qp, Inf)/norm(x_qp, Inf) <= 1e-5

println("Solving random SVM problem: fast proximal gradient (quiet)")
x, slv = solve(zeros(n), HingeLoss(b, mu), A, zeros(m), FPG(verbose = 0))
@test norm(x-x_qp, Inf)/norm(x_qp, Inf) <= 1e-5
