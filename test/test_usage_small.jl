A = randn(3,5)
b = randn(3)

x_pg = Variable(5)
prob_pg = problem(ls(A*x_pg - b) + norm(x_pg, 1))
sol_pg = solve!(prob_pg, PG())

x_zfpr = Variable(5)
prob_zfpr = problem(ls(A*x_zfpr - b) + norm(x_zfpr, 1))
sol_zfpr = solve!(prob_zfpr, ZeroFPR())
