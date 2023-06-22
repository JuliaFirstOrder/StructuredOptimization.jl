A = randn(3,5)
b = randn(3)

x_zfpr = Variable(5)
prob_zfpr = problem(ls(A*x_zfpr - b) + 1e-3*norm(x_zfpr, 1))
sol_zfpr = solve(prob_zfpr, ZeroFPR())

x_pnc = Variable(5)
prob_pnc = problem(ls(A*x_pnc - b) + 1e-3*norm(x_pnc, 1))
sol_pnc = solve(prob_pnc, PANOC())

x_pncp = Variable(5)
prob_pncp = problem(ls(A*x_pncp - b) + 1e-3*norm(x_pncp, 1))
sol_pncp = solve(prob_pncp, PANOCplus())
