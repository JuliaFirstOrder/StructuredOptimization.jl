module LineSpectraEstimation

using PyPlot
using RegLS
using Convex
using Mosek
using AbstractOperators
using DSP

function set_up()

	srand(17)

	Fs = 16e3
	Nt = 2^8 #time samples
	y0 = zeros(Nt)
	SNR = 10

	s = 6 #super-resolution factor
	fs = linspace(0,Fs/2,div(s*Nt,2)+1)            # super resolution frequency axis (up to Nyquist)
	t  = 0:1/Fs:(Nt-1)/Fs                          # time axis
	f  = linspace(0,Fs,Nt+1)[1:end-1]              # frequency axis
	K = 14                                         # number of sinusoids
	fk = fs[randperm(div(s*Nt,2)+1)[1:K]]          #sinusoids frequencies
	ak = 0.1*randn(K)+0.7                          # amplitude

	for i in eachindex(fk) y0 .+= ak[i].*sin.(2*π*fk[i].*t) end
	y = y0.+10^(-SNR/10)*sqrt(var(y0))*randn(length(y0))

	xzp = rfft([y;zeros((s-1)*length(y))])
	IDFTm = [exp(im*2*pi*k*n/(s*Nt))  for k =0:s*Nt-1, n=0:s*Nt-1] #Inverse Fourier Matrix
	S = [speye(Nt) spzeros(Nt,(s-1)*Nt)] # selection matrix
	Fm = S*IDFTm
	lambda_max_m = norm(Fm'*y, Inf)
	lambda_m = 0.06*lambda_max_m

	F = (s*Nt*IRDFT((div(s*Nt,2)+1,),s*Nt))[1:Nt] # Abstract Operator 
	lambda_max = norm(F'*y, Inf)
	lambda = 0.06*lambda_max

	setup = K, F, Fm, lambda, lambda_m
	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, setup

end

function solve_problem!(slv, x0, y, K, F, Fm, lambda, lambda_m)
	@minimize ls(F*x0-y)+lambda*norm(x0,1) with slv
	return x0
end

function solve_problem_ncvx!(slv, x0, y, K, F, Fm, lambda, lambda_m)
	@minimize ls(F*x0-y) st norm(x0,0) <= K with slv
	return x0
end

function solve_problem_matrix!(slv, x0, y, K, F, Fm, lambda, lambda_m)
	@minimize ls(Fm*x0-y)+lambda_m*norm(x0,1) with slv
	return x0
end

function solve_problem_ncvx_matrix!(slv, x0, y, K, F, Fm, lambda, lambda_m)
	@minimize ls(Fm*x0-y) st norm(x0,0) <= 2*K with slv
	return x0
end

function solve_problem_Convex!(slv, x0, y, K, F, Fm, lambda, lambda_m)
	problem = minimize(0.5*norm(Fm*x0-y,2)^2+lambda_m*norm(x0,1)) 
	return problem
end

function run_demo()

	t, f, fs, fk, ak, s, Nt, Fs, xzp, y, setup = set_up()

	x = RegLS.Variable(zeros(Complex{Float64},div(s*Nt,2)+1))
	slv = ZeroFPR()

	println("Solving LASSO (Abstract Operator)")
	@time solve_problem!(slv, x, y, setup...)
	x1 = copy(~x)

	println("Rafine solution by solving non-convex problem (Abstract Operator)")
	@time solve_problem_ncvx!(slv, x, y, setup...)
	x0 = copy(~x)

	x0m = RegLS.Variable(zeros(Complex{Float64},s*Nt))
	slv = ZeroFPR()

	println("Solving LASSO (Matrix Operator)")
	@time solve_problem_matrix!(slv, x0m, y, setup...)

	println("Rafine solution by solving non-convex problem (Matrix Operator)")
	@time solve_problem_ncvx_matrix!(slv, x0m, y, setup...)

	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, x0
end

function run_demo_Convex()
	t, f, fs, fk, ak, s, Nt, Fs, xzp, y, setup = set_up()
	x0m = Convex.ComplexVariable(s*Nt)

	println("Solving LASSO with Convex.jl (Matrix Operator)")
	slv = MosekSolver()
	problem = solve_problem_Convex!(slv, x0m, y, setup...)
	@time Convex.solve!(problem,slv)
	x1 = x0m.value
	x1 = x1[1:div(s*Nt,2)+1]
	return t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, zeros(x1)
end

function show_results(t, f, fs, fk, ak, s, Nt, Fs, xzp, y, x1, x0)

	figure()
	plot(fs,abs.(xzp./Nt ), label = "dft zero pad.")
	plot(f,        abs.(fft(y)./Nt)       , label = "dft")
	plot(fk,       abs.(ak)/2     , "r*", label = "true amp.")
	plot(fs,abs.(x1), "k*", label = "LASSO")
	plot(fs,abs.(x0), "go", label = "IndBallL0")
	xlim([0;Fs/2])
	legend()

end

end
