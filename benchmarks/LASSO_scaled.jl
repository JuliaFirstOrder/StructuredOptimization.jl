using StructuredOptimization
using BenchmarkTools
using JuMP
using SCS
using ECOS

function set_up(S::Int) #S scales problem

	n = S 
	m = div(S,5)
	SNR = 5

	srand(123)

	A = sprandn(m,n,5/n)
	println("n = $n, m = $m, nnz(A) = $(countnz(A))")
	x0 = zeros(n)
	x0[randperm(n)[1:div(n,4)+1]] = randn(div(n,4)+1)

	y = A*x0
	y += 10^(-SNR/10)*sqrt(var(y))*randn(length(y))
	lambda = 0.01*norm(A'*y,Inf) 
	x = StructuredOptimization.Variable(n)
	@minimize ls(A*x-y)+lambda*norm(x,1) with ZeroFPR(verbose = 0, tol =1e-12) 
	setup = A, y, lambda, n
	setup_JuMP = create_JuMP_model(A, y, lambda, n)
	return ~x, setup, setup_JuMP
end


function solve_problem(slv::S, A, y, lambda, n, M, xJ) where {S <: StructuredOptimization.ForwardBackwardSolver}
	x = StructuredOptimization.Variable(n)
	it, = @minimize ls(A*x-y)+lambda*norm(x,1) with slv
	return ~x, it, :UserLimit
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
	status = JuMP.solve(M)
	return getvalue(xJ), 0, status
end

function benchmark_LASSO(slv, eps_max, nvar)

	slv = eval(parse(slv))

	verbose, samples, seconds = 0, 5, 2*60
	xopt, setup, setup_JuMP = set_up(nvar)

	if slv == PG
		if nvars == 100000 #100000
			iterations = 1:div(100000,6):100000
		elseif nvars == 10000 #10000
			iterations = 1:div(50000,6):50000
		elseif nvars == 1000 #7875
			iterations = 1:div(7857,6):7857
		end
	elseif slv == FPG
		if nvars == 100000 #28941
			iterations = 1:div(28941,6):28941
		elseif nvars == 10000 #12534
			iterations = 1:div(13000,6):13000
		elseif nvars == 1000 #4027
			iterations = 1:div(4027,6):4027
		end
	elseif slv == ZeroFPR
		if nvars == 100000 #1019
			iterations = 1:div(1019,6):1019
		elseif nvars == 10000 #640
			iterations = 1:div(640,6):640
		elseif nvars == 1000 #163
			iterations = 1:div(163,6):163
		end
	elseif slv == PANOC
		if nvars == 100000 #1019
			iterations = 1:div(2*1019,6):2*1019
		elseif nvars == 10000 #640
			iterations = 1:div(2*640,6):2*640
		elseif nvars == 1000 #163
			iterations = 1:div(2*163,6):2*163
		end
	elseif slv == ECOSSolver
		if nvars == 100000 #104
			iterations = 1:div(104,6):104
		elseif nvars == 10000 #30
			iterations = 1:div(35,6):35
		elseif nvars == 1000 #23
			iterations = 1:div(25,6):25
		end
	elseif slv == SCSSolver
		if nvars == 100000 
			itSCS = 2680 #with default (low) tolerance
			iterations = 1:div(itSCS,3):6*div(itSCS,3)
		elseif nvars == 10000 
			itSCS = 420  #with default (low) tolerance
			iterations = 1:div(itSCS,3):8*div(itSCS,3)
		elseif nvars == 1000 
			itSCS = 180  #with default (low) tolerance
			iterations = 1:div(itSCS,3):6*div(itSCS,3)
		end
	end
		
	#iterations = [100000]
	err = zeros(length(iterations))
	t   = zeros(length(iterations))
	its = zeros(length(iterations))
	status = [:NotSolved]

	for i in eachindex(iterations) 

		if slv == SCSSolver
			solver = slv(verbose = verbose, max_iters = iterations[i], eps = 1e-20)
		elseif slv == ECOSSolver
			solver = slv(verbose = verbose, maxit = iterations[i])
		else
			solver = slv(verbose = verbose, maxit = iterations[i], tol = 1e-20)
		end

		b = @benchmarkable((x,it,status) = solve_problem(solver, setup..., setup_JuMP...), 
				   setup = (
					setup  = deepcopy($setup); 
					setup_JuMP  = deepcopy($setup_JuMP); 
					solver = deepcopy($solver);
					x = nothing;
					status = nothing;
					it = 0;
					), 
				   teardown = (
					   $err[$i] = norm(x-$xopt)/norm($xopt);
					   $status[1] = status;
					   ), 
				   evals = 1, samples = samples, seconds = seconds)

		results = run(b)
		t[i] = time(median(results))
		its[i] = iterations[i]
		if err[i] <= eps_max || status[1] == :Optimal
			t = t[1:i]
			err = err[1:i]
			its = its[1:i]
			break
		end
	end

	return t, err, its
end
BLAS.set_num_threads(5)

nvars = 100000
solvers = [
#	   "SCSSolver",
#  	   "ECOSSolver",
	   "PG",
	   "FPG",
	   "ZeroFPR",
	   "PANOC"
	   ]
T, ERR = [],[]
for slv in solvers
	t, err = benchmark_LASSO(slv, 1e-8, nvars)
	push!(T,t)
	push!(ERR,err)
end

using PyPlot
using DataFrames, CSV
save_stuff = true
dirpath = "/home/nantonel/Proximal_Gradient_Algorithms/fig"
#mkpath("data/")
figure()
for i in eachindex(solvers)
	plot(1e-9.*T[i], log10.(ERR[i]),  label = solvers[i],":*")
    if save_stuff
        CSV.write("$dirpath/lasso_$(nvars)_$(solvers[i]).cvs", 
                  DataFrame(q = log10.(ERR[i]), t = 1e-9.*T[i]) )
    end
end
legend()
ylabel("time")
xlabel("solution quality")


















