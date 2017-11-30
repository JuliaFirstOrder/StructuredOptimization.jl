
module SparseDeconvolution

using RegLS
using Convex
using Mosek
using AbstractOperators
using RIM
using PyPlot

function set_up()

	srand(123)
	Fs = 4000        # sampling frequency
	SNR = 15
	env = AcEnv(Fs)

	Lx,Ly,Lz = 4.,5.,3. # room dimensions
	T60 = 0.3           # reverberation time
	geo = CuboidRoom(Lx,Ly,Lz,T60,env) 

	xs = [0.5 0.5 0.5]'            # src pos (in meters)
	xr = [Lx-0.1 Ly-0.3 Lz-0.2]'   # mic pos
	Nh = div(Fs,4)                 # IR samples 
	Nx = div(Fs,2)                 # x  samples 

	t = linspace(0,1/Fs*Nx,Nx)

	h = rim(xs,xr,Nh,geo,env)[:] # Impuse Response

	x = full(sprandn(Nx, 0.05))  # sparse input 

	y0 = conv(x,h)                                      # output signal
	y = y0+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0)) # add noise

	H = Conv(Float64,size(x),h)                              # Abstract Operator
	T = hcat([[zeros(i);h;zeros(Nx-1-i)] for i = 0:Nx-1]...) # Full Matrix
	lambda = 1e-2*vecnorm(H'*y,Inf)

	setup = y, H, T, lambda 
	return setup, t, x

end

function solve_problem!(slv, x0, y, H, T, lambda)
	@minimize ls(H*x0-y)+lambda*norm(x0,1) with slv
	return x0
end

function solve_problem_matrix!(slv, x0, y, H, T, lambda)
	@minimize ls(T*x0-y)+lambda*norm(x0,1) with slv
	return x0
end

function solve_problem_Convex!(slv, x0, y, H, T, lambda)
	problem = minimize(0.5*norm(T*x0-y,2)^2+lambda*norm(x0,1)) 
	return problem
end


function run_demo()

	setup, t, x = set_up()
	x0 = RegLS.Variable(zeros(x))

	println("Solving Unregularized problem")
	xu = setup[3]\setup[1]

	println("Solving Regularized problem with Abstract Operator")
	slv = ZeroFPR()
	@time x0 = solve_problem!(slv, x0, setup...)
	x1 = copy(~x0)

	println("Solving Regularized problem with Full Matrix")
	~x0 .= 0
	slv = ZeroFPR()
	@time x0 = solve_problem_matrix!(slv, x0, setup...)
	xm = copy(~x0)

	println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
	println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

	return t, x, x1, xu
end

function run_demo_Convex()

	setup, t, x = set_up()
	x0 = Convex.Variable(size(x)...)

	println("Solving Unregularized problem")
	xu = setup[3]\setup[1]

	println("Solving Regularized problem with Convex")
	slv = MosekSolver()
	problem = solve_problem_Convex!(slv, x0, setup...)
	@time Convex.solve!(problem,slv)
	x1 = x0.value

	println("  regularized MSE: $( 20*log10(norm(x1-x)/norm(x)) )")
	println("unregularized MSE: $( 20*log10(norm(xu-x)/norm(x)) )")

	return t, x, x1, xu
end

function show_results(t, x, x1, xu)

	figure()
	plot(t,x )
	plot(t,x1)
	plot(t,xu)

end

end

