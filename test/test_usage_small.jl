A = randn(3,5)
b = randn(3)

x_pg = Variable(5)
prob_pg = problem(ls(A*x_pg - b) + 1e-3*norm(x_pg, 1))
sol_pg = solve(prob_pg, StructuredOptimization.PG())

x_zfpr = Variable(5)
prob_zfpr = problem(ls(A*x_zfpr - b) + 1e-3*norm(x_zfpr, 1))
sol_zfpr = solve(prob_zfpr, StructuredOptimization.ZeroFPR())

x_pnc = Variable(5)
prob_pnc = problem(ls(A*x_pnc - b) + 1e-3*norm(x_pnc, 1))
sol_pnc = solve(prob_pnc, StructuredOptimization.PANOC())
