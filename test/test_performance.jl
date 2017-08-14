A = randn(100, 200)
B = randn(100, 300)
c = randn(100)

x_pg = Variable(200)
y_pg = Variable(300)
prob_pg = problem(ls(A*x_pg + B*y_pg - c) + norm(x_pg, 1) + norm(y_pg))
sol_pg = solve!(prob_pg, PG(verbose=0))

x_pg = Variable(200)
y_pg = Variable(300)
@time prob_pg = problem(ls(A*x_pg + B*y_pg - c) + norm(x_pg, 1) + norm(y_pg))
@time sol_pg = solve!(prob_pg, PG(verbose=0))

println(sol_pg)

x_fpg = Variable(200)
y_fpg = Variable(300)
prob_fpg = problem(ls(A*x_fpg + B*y_fpg - c) + norm(x_fpg, 1) + norm(y_fpg))
sol_fpg = solve!(prob_fpg, FPG(verbose=0))

x_fpg = Variable(200)
y_fpg = Variable(300)
@time prob_fpg = problem(ls(A*x_fpg + B*y_fpg - c) + norm(x_fpg, 1) + norm(y_fpg))
@time sol_fpg = solve!(prob_fpg, FPG(verbose=0))

println(sol_fpg)

x_zfpr = Variable(200)
y_zfpr = Variable(300)
prob_zfpr = problem(ls(A*x_zfpr + B*y_zfpr - c) + norm(x_zfpr, 1) + norm(y_zfpr))
sol_zfpr = solve!(prob_zfpr, ZeroFPR(verbose=0))

x_zfpr = Variable(200)
y_zfpr = Variable(300)
@time prob_zfpr = problem(ls(A*x_zfpr + B*y_zfpr - c) + norm(x_zfpr, 1) + norm(y_zfpr))
@time sol_zfpr = solve!(prob_zfpr, ZeroFPR(verbose=0))

println(sol_zfpr)
