module MatrixDecomposition

using BenchmarkTools
using StructuredOptimization
using Images, ImageView

function set_up()
    # dataset from: http://pione.dinf.usherbrooke.ca/static/dataset/clutter/IndianTraffic3.zip
	Frames = 900
	dir = "street"
	img = load("$dir/in000000.jpg")
	n,m = size(img,1),size(img,2)
    frames = collect(1:135:Frames)
	N = length(frames)

	Y = zeros(Float64,n,m,N)
	for f in eachindex(frames)
		a = @sprintf("%6.6i",frames[f])
		img = load("$dir/in$(a).jpg")
		Y[:,:,f] .= convert(Array{Float64},Gray.(img))
		img = []
	end

	Y = reshape(Y,n*m,N)

	F = Variable(n*m,N)
	B = Variable(n*m,N)

    println("frames: $N, number of variables: $(2*N*m*n)")

	R, lambda = 1, 3e-2
	return B, F, Y, R, lambda, n, m, N 
end

function run_demo()
	slv = PANOC(tol = 1e-4, verbose = 2)
	setup = set_up()
	solve_problem!(slv,setup...)
	return setup
end

function solve_problem!(slv, B, F, Y, R, lambda, n, m, N)
	it, = @minimize ls(B+F-Y)+lambda*norm(F,1) st rank(B) <= R with slv 
	return it
end

function benchmark(;verb = 0, samples = 5, seconds = 100, tol = 1e-4)

	suite = BenchmarkGroup()
	
	solvers = [
             #  "ZeroFPR",
               "PANOC",
               "PG"
              ]
	slv_opt = [
              # "(verbose = $verb, tol = $tol)", 
               "(verbose = $verb, tol = $tol)",
               "(verbose = $verb, tol = $tol)"
              ]

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
	println("MatrixDecomposition its")
	println(its)
	return results
end

function show_results(B, F, Y, R, lambda, n, m, N)
	F = F.x
	B = B.x

	F[F .!=0] = F[F.!=0]+reshape(B,size(F))[F.!=0]
	F[F.== 0] = 1. #put white in null pixels
	F[F.>1] .= 1; F[F.<0.] .= 0.
	B[B.>1] .= 1; B[B.<0.] .= 0.

	#convert back to images
	Y = reshape(Y,n,m,N)
	B = reshape(B,n,m,N)
	F = reshape(F,n,m,N)
	#
	####videos
	#imshow(Y)
	#imshow(F)
	#imshow(B)
	#
	##images
	idx = round.(Int,linspace(1,N,4)) 
	imshow(hcat([Y[:,:,i] for i in idx]...))
	imshow(hcat([F[:,:,i] for i in idx]...))
	imshow(hcat([B[:,:,i] for i in idx]...))

    save_stuff = false
    if save_stuff == true
        for (i,idxi) in enumerate(idx)
            save("Y$i.jpg", Y[:,:,idxi])
            save("F$i.jpg", F[:,:,idxi])
            save("B$i.jpg", B[:,:,idxi])
        end
    end
    
end


end
