@printf("\nTesting cost terms\n")

# Simple Terms

x = Variable(10)
cf = norm(x, 0)
cf = norm(x, 0) <= 3
cf = norm(x, 1)
cf = norm(x, 1) <= 1.5
cf = norm(x, 2)
cf = norm(x, 2) <= 2.3
cf = norm(x, Inf)
cf = x <= 3.0
cf = 3.0 <= x
cf = x >= 1.0
cf = 1.0 >= x
cf = x in [-5.0, 5.0]
cf = norm(x)^2
cf = norm(x, 2)^2

X = Variable(10, 10)
cf = rank(X) <= 3

# Summing terms

x = Variable(10)
cf = ls(x) + norm(x, 1)

# More complex situations

x = Variable(10)
A = randn(5, 10)
y = Variable(7)
B = randn(5, 7)
b = randn(5)

cf = ls(A*x - b) + norm(x, 1)
cf1 = ls(A*x - B*y + b) + norm(y, 1) + norm(y, 2)
cf2 = 0.5*norm(A*x - B*y + b, 2)^2 + norm(y, 1) + norm(y, 2)
