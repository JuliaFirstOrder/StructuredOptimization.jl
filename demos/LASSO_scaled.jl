
using RegLS
using BenchmarkTools
using JuMP
using SCS
using Mosek
using ECOS

function set_up(S::Int) #S scales problem

	n = S 
	m = div(S,4)
	SNR = 5

	srand(123)

	A = sprandn(m,n,5/n)
	println("n = $n, nnz(A) = $(countnz(A))")
	x0 = zeros(n)
	x0[randperm(n)[1:div(n,4)+1]] = randn(div(n,4)+1)
	x0 /= norm(x0)^2

	y = A*x0
	y += 10^(-SNR/10)*sqrt(var(y))*randn(length(y))
	lambda = 0.01*norm(A'*y,Inf) 
	x = RegLS.Variable(n)
	@minimize ls(A*x-y)+lambda*norm(x,1) with ZeroFPR(verbose = 0, tol =1e-12) 
	setup = A, y, lambda, n
	setup_JuMP = create_JuMP_model(A, y, lambda, n)
	return ~x, setup, setup_JuMP
end


function solve_problem(slv::S, A, y, lambda, n, M, xJ) where {S <: RegLS.ForwardBackwardSolver}
	x = RegLS.Variable(n)
	slv = @minimize ls(A*x-y)+lambda*norm(x,1) with slv
	return ~x, slv.it
end

function create_JuMP_model(A, y, lambda, n)
	M = Model()
	@variables M begin
		x[1:n]
		t[1:n]
		w
	end
	@objective(M,Min,[0.5;lambda*ones(n)]'*[w;t])
	@constraint(M, soc, norm( [1-w;2*(A*x-y)] ) <= 1+w)
	@constraint(M,  x .<= t)
	@constraint(M, -t .<= x)
	return M, x
end

function solve_problem(slv::S, A, y, lambda, n, M, xJ) where {S <: MathProgBase.SolverInterface.AbstractMathProgSolver}
	setsolver(M, slv)
	solve(M)
	return getvalue(xJ), 0
end

function benchmark_LASSO()
	suite = BenchmarkGroup()

	verbose, samples, seconds = 0, 5, 30*60

	solvers = [
		   "ECOSSolver",
		   "SCSSolver", 
		   "PG", 
		   "FPG", 
		   "ZeroFPR", 
		   ]
	slv_opt = [
		   "(verbose = $verbose, maxit     = 10000)", 
		   "(verbose = $verbose, max_iters = 100000)", 
		   "(verbose = $verbose, maxit     = 100000, tol = 1e-6)", 
		   "(verbose = $verbose, maxit     = 100000, tol = 1e-6)", 
		   "(verbose = $verbose, maxit     = 100000, tol = 1e-6)"
		   ]
	iterations = Dict([(sol,0) for sol in solvers]) 
	nvar =  [1000;10000;100000]
	#nvar =  [100]
	err = Dict((n,Dict([(sol,0.) for sol in solvers])) for n in nvar) 
	its = Dict((n,Dict([(sol,0.) for sol in solvers])) for n in nvar) 

	for ii in nvar 
		suite[ii] = BenchmarkGroup()
		xopt, setup, setup_JuMP = set_up(ii)
		for i in eachindex(solvers)
			solver = eval(parse(solvers[i]*slv_opt[i]))

			suite[ii][solvers[i]] = 
			@benchmarkable((x,it) = solve_problem(solver, setup..., setup_JuMP...), 
				       setup = (
						setup  = deepcopy($setup); 
						setup_JuMP  = deepcopy($setup_JuMP); 
						solver = deepcopy($solver);
						x = nothing;
						it = 0;
						), 
				       teardown = (
						   $err[$ii][$solvers[$i]] = norm(x-$xopt);
						   $its[$ii][$solvers[$i]] = it;
						   ), 
				       evals = 1, samples = samples, seconds = seconds)
		end

	end

	benchmarks = run(suite)

	return benchmarks, solvers, err, its, nvar
end
BLAS.set_num_threads(4)

benchmarks, solvers, err, its, nvar = benchmark_LASSO()

println("\n")
showall(median(benchmarks))
println("\n")

#using PyPlot
#figure()
#for slv in solvers
#	semilogx(nvar,[10*log10(time(median(benchmarks[i][slv]))) for i in nvar], label = slv)
#end
#legend()

import BenchmarkTools:prettytime
tab = "\\midrule \n"
for i in nvar
	tab *= "\\multirow{2}{*}{ \$ n = 10^$(Int(log10(i))) \$ } & "
	tab *= "\$ t \$ "
	for slv in solvers 
		tab *= "&  $(prettytime(time(median(benchmarks[i][slv]))))  "
	end
	tab *= "\\\\ \n\\cmidrule(lr){2-7}\n                                & \$ \\epsilon \$ " 
	for slv in solvers 
		tab *= "& $(round(20*log10(err[i][slv]),1)) "
	end
	tab *= "\\\\ \n \\midrule \n" 
end

tab = replace(tab, "μ", "\$\\mu\$")
println(tab)
write("table.tex", tab)

















