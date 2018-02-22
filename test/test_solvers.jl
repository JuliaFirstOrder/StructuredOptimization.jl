@printf("\nTesting solvers\n")

srand(0)
tol = 1e-6;

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

normL = norm(A1)^2

# Apply PG

x = zeros(n)
sol = StructuredOptimization.apply!(StructuredOptimization.PG(tol=tol), x; fq=f, Aq=L, g=g)

gstep = x - (A1'*(A1*x-b))
pgstep = sign.(gstep).*max.(0, abs.(gstep) .- lam)
@test norm(pgstep - x)/normL^2 <= tol
# project subgr onto the subdifferential of the L1-norm at x
subgr = -A1'*(A1*x-b)
subgr_proj = min.(max.(subgr, -lam), lam)
subgr_proj[x .< 0] = -lam
subgr_proj[x .> 0] = lam
@test norm(subgr - subgr_proj, Inf) <= 1e-5

# Apply FPG

x = zeros(n)
sol = StructuredOptimization.apply!(StructuredOptimization.FPG(tol=tol), x; fq=f, Aq=L, g=g)

gstep = x - (A1'*(A1*x-b))
pgstep = sign.(gstep).*max.(0, abs.(gstep) .- lam)
@test norm(pgstep - x)/normL^2 <= tol
# project subgr onto the subdifferential of the L1-norm at x
subgr = -A1'*(A1*x-b)
subgr_proj = min.(max.(subgr, -lam), lam)
subgr_proj[x .< 0] = -lam
subgr_proj[x .> 0] = lam
@test norm(subgr - subgr_proj, Inf) <= 1e-5

# Apply ZeroFPR

x = zeros(n)
zerofpr = StructuredOptimization.ZeroFPR(tol=tol,maxit=1000,verbose=1)
sol = StructuredOptimization.apply!(zerofpr, x; fq=f, Aq=L, g=g)

gstep = x - (A1'*(A1*x-b))
pgstep = sign.(gstep).*max.(0, abs.(gstep) .- lam)
@test norm(pgstep - x)/normL^2 <= tol
# project subgr onto the subdifferential of the L1-norm at x
subgr = -A1'*(A1*x-b)
subgr_proj = min.(max.(subgr, -lam), lam)
subgr_proj[x .< 0] = -lam
subgr_proj[x .> 0] = lam
@test norm(subgr - subgr_proj, Inf) <= 1e-5

# Apply PANOC

x = zeros(n)
panoc = StructuredOptimization.PANOC(tol=tol,maxit=1000,verbose=1)
sol = StructuredOptimization.apply!(panoc, x; fq=f, Aq=L, g=g)

gstep = x - (A1'*(A1*x-b))
pgstep = sign.(gstep).*max.(0, abs.(gstep) .- lam)
@test norm(pgstep - x)/normL^2 <= tol
# project subgr onto the subdifferential of the L1-norm at x
subgr = -A1'*(A1*x-b)
subgr_proj = min.(max.(subgr, -lam), lam)
subgr_proj[x .< 0] = -lam
subgr_proj[x .> 0] = lam
@test norm(subgr - subgr_proj, Inf) <= 1e-5

###########################################################################
# Regularized/constrained least squares with two variable blocks
###########################################################################

m, n, l = 10, 30, 40
A1 = randn(m, n)
A2 = randn(m, l)
x1 = randn(n)
x2 = randn(l)
x2 = 1.1.*x2./norm(x2)
b = A1*x1 + A2*x2

f = PrecomposeDiagonal(SqrNormL2(), 1.0, -b)
g = SeparableSum((NormL1(lam), IndBallL2(1.0)))
L = HCAT(MatrixOp(A1), MatrixOp(A2))

normL = norm([A1 A2])

# Apply PG

x = (zeros(n), zeros(l))
sol = StructuredOptimization.apply!(StructuredOptimization.PG(tol=tol), x; fq=f, Aq=L, g=g)

res = A1*x[1]+A2*x[2]-b
gstep1 = x[1] - A1'*res
gstep2 = x[2] - A2'*res
pgstep1 = sign.(gstep1).*max.(0, abs.(gstep1) .- lam)
pgstep2 = norm(gstep2) > 1 ? gstep2/norm(gstep2) : gstep2
@test norm(x[1] - pgstep1)/normL^2 <= tol
@test norm(x[2] - pgstep2)/normL^2 <= tol

# Apply FPG

x = (zeros(n), zeros(l))
sol = StructuredOptimization.apply!(StructuredOptimization.FPG(tol=tol), x; fq=f, Aq=L, g=g)

res = A1*x[1]+A2*x[2]-b
gstep1 = x[1] - A1'*res
gstep2 = x[2] - A2'*res
pgstep1 = sign.(gstep1).*max.(0, abs.(gstep1) .- lam)
pgstep2 = norm(gstep2) > 1 ? gstep2/norm(gstep2) : gstep2
@test norm(x[1] - pgstep1)/normL^2 <= tol
@test norm(x[2] - pgstep2)/normL^2 <= tol

# Apply ZeroFPR

x = (zeros(n), zeros(l))
sol = StructuredOptimization.apply!(StructuredOptimization.ZeroFPR(tol=tol), x; fq=f, Aq=L, g=g)

res = A1*x[1]+A2*x[2]-b
gstep1 = x[1] - A1'*res
gstep2 = x[2] - A2'*res
pgstep1 = sign.(gstep1).*max.(0, abs.(gstep1) .- lam)
pgstep2 = norm(gstep2) > 1 ? gstep2/norm(gstep2) : gstep2
@test norm(x[1] - pgstep1)/normL^2 <= tol
@test norm(x[2] - pgstep2)/normL^2 <= tol

# Apply PANOC

x = (zeros(n), zeros(l))
sol = StructuredOptimization.apply!(StructuredOptimization.PANOC(tol=tol), x; fq=f, Aq=L, g=g)

res = A1*x[1]+A2*x[2]-b
gstep1 = x[1] - A1'*res
gstep2 = x[2] - A2'*res
pgstep1 = sign.(gstep1).*max.(0, abs.(gstep1) .- lam)
pgstep2 = norm(gstep2) > 1 ? gstep2/norm(gstep2) : gstep2
@test norm(x[1] - pgstep1)/normL^2 <= tol
@test norm(x[2] - pgstep2)/normL^2 <= tol

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
sol = StructuredOptimization.apply!(StructuredOptimization.PG(tol=tol), x; fq=f, Aq=L, g=g)

# Apply FPG

x = zeros(n)
sol = StructuredOptimization.apply!(StructuredOptimization.FPG(tol=tol), x; fq=f, Aq=L, g=g)

# Apply ZeroFPR

x = zeros(n)
sol = StructuredOptimization.apply!(StructuredOptimization.ZeroFPR(tol=tol), x; fq=f, Aq=L, g=g)

# Apply PANOC

x = zeros(n)
sol = StructuredOptimization.apply!(StructuredOptimization.PANOC(tol=tol), x; fq=f, Aq=L, g=g)
