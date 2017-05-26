@printf("\nTesting solvers\n")

###########################################################################
# Lasso
###########################################################################

m, n, l = 10, 50, 30
A1 = randn(m, n)
A2 = randn(m, l)
b = randn(m)
lam_max = norm(A1'*b, Inf)
lam = 0.05*lam_max

f = PrecomposeDiagonal(SqrNormL2(), 1.0, -b)
g = NormL1(lam)
L = MatrixOp(A1)

# Apply PG

x = zeros(n)
sol = RegLS.apply(PG(), x, f, L, g)

gstep = x - (A1'*(A1*x-b))
pgstep = sign.(gstep).*max.(0, abs.(gstep) .- lam)
@test norm(pgstep - x, Inf) <= 1e-6
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
@test norm(pgstep - x, Inf) <= 1e-6
# project subgr onto the subdifferential of the L1-norm at x
subgr = -A1'*(A1*x-b)
subgr_proj = min.(max.(subgr, -lam), lam)
subgr_proj[x .< 0] = -lam
subgr_proj[x .> 0] = lam
@test norm(subgr - subgr_proj, Inf) <= 1e-5
