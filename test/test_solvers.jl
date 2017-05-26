@printf("\nTesting solvers\n")

###########################################################################
# Lasso
###########################################################################

m, n, l = 10, 50, 30
A1 = randn(m, n)
b = randn(m)
lam_max = norm(A1'*b, Inf)
lam = 0.3*lam_max

f = PrecomposeDiagonal(SqrNormL2(), 1.0, -b)
g = NormL1(lam)
L = MatrixOp(A1)

# Apply PG

x = zeros(n)
sol = RegLS.apply(PG(), x, f, L, g)

gstep = x - (A1'*(A1*x-b))
pgstep = sign.(gstep).*max.(0, abs.(gstep) .- lam)
@test norm(pgstep - x) <= 1e-8
# project subgr onto the subdifferential of the L1-norm at x
subgr = -A1'*(A1*x-b)
subgr_proj = min.(max.(subgr, -lam), lam)
subgr_proj[x .< 0] = -lam
subgr_proj[x .> 0] = lam
@test norm(subgr - subgr_proj, Inf) <= 1e-5

# Apply FPG

x = zeros(n)
sol = RegLS.apply(FPG(), x, f, L, g)

gstep = x - (A1'*(A1*x-b))
pgstep = sign.(gstep).*max.(0, abs.(gstep) .- lam)
@test norm(pgstep - x) <= 1e-8
# project subgr onto the subdifferential of the L1-norm at x
subgr = -A1'*(A1*x-b)
subgr_proj = min.(max.(subgr, -lam), lam)
subgr_proj[x .< 0] = -lam
subgr_proj[x .> 0] = lam
@test norm(subgr - subgr_proj, Inf) <= 1e-5

###########################################################################
# Regularized/constrained least squares with two variable blocks
###########################################################################

m, n, l = 10, 50, 30
A1 = randn(m, n)
A2 = randn(m, l)
x1 = randn(n)
x2 = randn(l)
x2 = 1.2.*x2./norm(x2)
b = A1*x1 + A2*x2

f = PrecomposeDiagonal(SqrNormL2(), 1.0, -b)
g = SeparableSum((NormL1(lam), IndBallL2(1.0)))
L = HCAT(MatrixOp(A1), MatrixOp(A2))

# Apply PG

x = (zeros(n), zeros(l))
sol = RegLS.apply(PG(), x, f, L, g)

res = A1*x[1]+A2*x[2]-b
gstep1 = x[1] - A1'*res
gstep2 = x[2] - A2'*res
pgstep1 = sign.(gstep1).*max.(0, abs.(gstep1) .- lam)
pgstep2 = norm(gstep2) > 1 ? gstep2/norm(gstep2) : gstep2
@test norm(x[1] - pgstep1) <= 1e-8
@test norm(x[2] - pgstep2) <= 1e-8

# Apply FPG

x = (zeros(n), zeros(l))
sol = RegLS.apply(FPG(), x, f, L, g)

res = A1*x[1]+A2*x[2]-b
gstep1 = x[1] - A1'*res
gstep2 = x[2] - A2'*res
pgstep1 = sign.(gstep1).*max.(0, abs.(gstep1) .- lam)
pgstep2 = norm(gstep2) > 1 ? gstep2/norm(gstep2) : gstep2
@test norm(x[1] - pgstep1) <= 1e-8
@test norm(x[2] - pgstep2) <= 1e-8

###########################################################################
# L2-regularized least squares with two data blocks (just to play)
###########################################################################

m1, m2, n, = 30, 40, 50
A1 = randn(m1, n)
A2 = randn(m2, n)
b1 = randn(m1)
b2 = randn(m2)

f1 = PrecomposeDiagonal(SqrNormL2(), 1.0, -b1)
f2 = PrecomposeDiagonal(SqrNormL2(), 1.0, -b2)
f = SeparableSum((f1, f2))
g = NormL2(1.0)
L = VCAT(MatrixOp(A1), MatrixOp(A2))

# Apply PG

x = zeros(n)
sol = RegLS.apply(PG(), x, f, L, g)

# Apply FPG

x = zeros(n)
sol = RegLS.apply(FPG(), x, f, L, g)
