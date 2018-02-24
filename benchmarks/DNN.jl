module DNN

using BenchmarkTools
using StructuredOptimization, AbstractOperators, ProximalOperators
using PyPlot

function dnn(D::Matrix, n::Int, lambda)

	dim_in, N = size(D)
	W1 = Variable(sqrt(2/(2+n)*N)*randn(n,dim_in)) # first layer
	S1 = Sigmoid((n,N)) 
	b1 = Variable(1)
	W2 = Variable(sqrt(2/(2*n))*randn(n,n))        # second layer
	S2 = Sigmoid((n,N)) 
	b2 = Variable(1)
	W3 = Variable(sqrt(2/(n+1))*randn(1,n))        # third layer
	S3 = Sigmoid((1,N)) 
	b3 = Variable(1)

	L1 = S1*(W1* D.+b1)
	L2 = S2*(W2*L1.+b2)
	y  = S3*(W3*L2.+b3)

	lambda1 = lambda/(2*n)
	lambda2 = lambda/(n+n)
	lambda3 = lambda/n

	reg = lambda1*norm(W1,1)+lambda2*norm(W2,1)+lambda3*norm(W3,1)
	return y, reg 
end

function create_dataset()

	srand(11)
	#construct dataset
	Na, Nb = 800,800
	N = Na+Nb  
	rhoa = 5.
	a = hcat([ rhoa.*[ cos(theta);  sin(theta)] for theta in linspace(0,2*pi,Na) ]...)   
	a .+= 1.5.*randn(size(a))
	b = zeros(2, div(Nb,2))
	rhob = 8.
	b = [b hcat([ rhob.*[ cos(theta);  sin(theta)] for theta in linspace(0,2*pi,div(Nb,2)) ]...)]   
	b .+= randn(size(b))
	yt = [ones(Bool,Na); zeros(Bool,Nb)] #labels

	A = [a b]
	#preprocess data
	A .-= mean(A)
	A ./= sqrt.(var(A,2))

	return Na, Nb, A, yt

end

function set_up()
	Na, Nb, A, yt = create_dataset()
	KCV = false
	# training
	K = 5  # KCV folds
	L = 5   # number of candidate lambdas
	#
	n = 7  # inner layers nodes
	if KCV == true  #K-fold cross validation
		lambda = cross_validation(Na, Nb, A, yt, n, K, L)
	else
		lambda = 0.0031622776601683794 # set lambda manually/already trained 
	end
	#create deep neural network
	y, reg = dnn(A, n, lambda)
	return Na, n, y, yt, A, reg
end

function run_demo()
	slv = PANOC(tol = 1e-4, maxit = 50000)
	setup = set_up()
	solve_problem!(slv,setup...)
	return setup
end

function cross_validation(Na,Nb,A,yt, n,K,L)


	Nacv, Nbcv = div(Na,K), div(Nb,K)
	fa,fb = reshape(randperm(K*Nacv),Nacv,K),reshape(randperm(K*Nbcv),Nbcv,K) #fold indices

	ytcv = [ones(Bool,Nacv);zeros(Bool,Nbcv)]              #labels validation
	yttr = [ones(Bool,(K-1)*Nacv);zeros(Bool,(K-1)*Nbcv)]  #labels training

	#initialize data matrices
	Dcv = zeros(2,Nacv+Nbcv)
	Dtr = zeros(2,(K-1)*(Nacv+Nbcv))
		
	ycv,     = dnn(Dcv,n,1.)
	ytr, reg = dnn(Dtr,n,1.)
	 
	lambdas = logspace(-4,-1,L)

	cv_err = zeros(K,L) # cross validation error
	ft_err = zeros(K,L) # fit error
	@views a, b = A[:,1:Na], A[:,Na+1:end]

	for k = 1:K, l = 1:L 

		if l == 1
			@printf("|------------------------------------------------------------------------| \n")
			ytr, reg = dnn(Dtr,n,1.) #re-initialize weights
		end

		@views Dcv .= [a[:,fa[:,k]] b[:,fb[:,k]]]                 # validation set
		idx_tr = find(1:K .!= k)
		@views Dtr .= [a[:,fa[:,idx_tr][:]] b[:,fb[:,idx_tr][:]]] # training set

		slv = PANOC(verbose = 0, tol = 1e-4, maxit = 5000) 
		@minimize crossentropy(ytr,yttr)+lambdas[l]*reg with slv 

		ft_err[k,l] = CrossEntropy(yttr)(operator(ytr)*(~).(variables(ytr)))
		cv_err[k,l] = CrossEntropy(ytcv)(operator(ycv)*(~).(variables(ytr)))

		@printf("| fold %2.1d / %2.1d | λ = %4.4f | %2.1d / %2.1d | ϵ cv =  %4.4f  | ϵ ft =  %4.4f | \n", 
			k, K, lambdas[l], l, L, cv_err[k,l], ft_err[k,l] )
	end

	println("trained lambda: $( lambdas[indmin(mean(cv_err,1))] )")
	#use lambda that gives minimum cross validation error
	return lambdas[indmin(mean(cv_err,1))]
end

		
function solve_problem!(slv, Na, n, y, yt, A, reg)
	it, = @minimize crossentropy(y,yt)+reg with slv 
	return it
end

function benchmark(;verb = 0, samples = 5, seconds = 100, tol = 1e-4, maxit = 50000 )

	suite = BenchmarkGroup()

	solvers = ["ZeroFPR", 
               "PANOC",
               "PG"]
	slv_opt = [
               "(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)", 
               "(verbose = $verb, tol = $tol, maxit = $maxit)"]

	its = Dict([(sol,0.) for sol in solvers])
	for i in eachindex(solvers)

		setup = set_up()
		solver = eval(parse(solvers[i]*slv_opt[i]))

		suite[solvers[i]] = 
		@benchmarkable(it = solve_problem!(solver, setup...), 
			       setup = ( 
					it = 0;
					setup = deepcopy($setup); 
					solver = deepcopy($solver) ), 
			       teardown = (
					  $its[$solvers[$i]] = it;
					  ), 
			       evals = 1, samples = samples, seconds = seconds)
	end

	results = run(suite, verbose = (verb != 0))
	println("DNN its")
	println(its)
	return results
end


function show_results(Na, n, y, yt, A, reg)

	xx = linspace(-3,3,200)
	yy = linspace(-3,3,200)

	R = vcat([[xi yi] for yi in yy, xi in xx][:]...)'
	yr,        = dnn(R,n,1.)
	regions = operator(yr)*(~).(variables(y))
	regions = reshape(regions,200,200)

	figure()
	pcolormesh(xx,yy,regions) 
	cnt = contour(xx,yy,regions, levels = 0.7*[maximum(regions)] ) 
	plot(A[1,Na+1:end],A[2,Na+1:end], "ko")
	plot(A[1,1:Na],A[2,1:Na], "r*")
end

end
