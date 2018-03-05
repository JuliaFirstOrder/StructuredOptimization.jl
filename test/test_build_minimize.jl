@printf("\n Testing solver build \n")

x = Variable(10)
A = randn(5, 10)
y = Variable(7)
B = randn(5, 7)
b = randn(5)

prob = problem(ls(A*x + b), norm(x, 2) <= 1.0)
built_slv = build(prob, StructuredOptimization.PG())
solve!(built_slv)

~x .= 0.
prob = problem(ls(A*x - B*y + b) + norm(y, 1), norm(x, 2) <= 1.0)
built_slv = build(prob, FPG())
solve!(built_slv)

@printf("\n Testing @minimize \n")
~x .= 0.
~y .= 0.
slv, = @minimize ls(A*x - B*y + b) st norm(x, 2) <= 1e4, norm(y, 1) <= 1.0 with PG()
~x .= 0.
slv, = @minimize ls(A*x - b) st norm(x, 1) <= 1.0 with PG()
~x .= 0.
slv, = @minimize ls(A*x - b) st norm(x, 1) <= 1.0
~x .= 0.
slv, = @minimize ls(A*x - b) + norm(x, 1) with PG()
~x .= 0.
slv, = @minimize ls(A*x - b) + norm(x, 1)
~x .= 0.
slv, = @minimize ls(A*x - b)

#TODO many many more tests
x = Variable(5)
A = randn(10, 5)
b = randn(10)

@printf("\n Testing @minimize nonlinear \n")
slv, = @minimize ls(sigmoid(A*x,10) - b)+norm(x,1) with PG()
xpg = copy(~x)
~x .= 0.
slv, = @minimize ls(sigmoid(A*x,10) - b)+norm(x,1) with ZeroFPR()
xz = copy(~x)
~x .= 0.
slv, = @minimize ls(sigmoid(A*x,10) - b)+norm(x,1) with PANOC()
xp = copy(~x)
~x .= 0.

@test norm(xz-xpg) <1e-4
@test norm(xp-xpg) <1e-4
