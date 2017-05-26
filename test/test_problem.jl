@printf("\nTesting problem construction\n")

x = Variable(10)
A = randn(5, 10)
y = Variable(7)
B = randn(5, 7)
b = randn(5)

prob = problem(0.5*norm(A*x - B*y + b, 2)^2 + norm(y, 1), norm(x, 2) <= 1.0)

sol = solve(prob, PG())
