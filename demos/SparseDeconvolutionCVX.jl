
#same using Convex
module SparseDeconvolutionCVX

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
	Nh = div(Fs,5)                 # IR samples 
	Nx = div(Fs,5)                 # x  samples 

	t = linspace(0,1/Fs*Nx,Nx)

	h = rim(xs,xr,Nh,geo,env)[:] # Impuse Response

	x = full(sprandn(Nx, 0.05))  # sparse input 

	y0 = conv(x,h)                                      # output signal
	y = y0+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0)) # add noise

	H = Conv(Float64,size(x),h)                              # Abstract Operator
	T = hcat([[zeros(i);h;zeros(Nx-1-i)] for i = 0:Nx-1]...) # Full Matrix
	x0 = Variable(size(x)...)
	lambda = 1e-2*vecnorm(H'*y,Inf)

	problem = minimize(sumsquares(T*x0-y)+lambda*norm_1(x0))

	return t, y, x, x0, problem

end

function solve_problem!(problems, solver)
	solve!(problems[1],solver)
end


function run_demo()
	t, y, x, x0,  = set_up()

	println("Solving Unregularized problem")
	xu = setup[3]\setup[1]

	println("Solving Regularized problem with Convex")
	x0.value = zeros(x)
	solver = MosekSolver()
	@time solve_problem!(setup)
	x1 = copy(x0.value)

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
