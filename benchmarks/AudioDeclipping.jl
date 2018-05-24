module AudioDeclipping

using BenchmarkTools
using StructuredOptimization, AbstractOperators
using DSP, WAV
using PyPlot

function set_up()

	Nl = 2^10                           # frames length
	Fs = 16e3                           # sampling frequency
	x0, Fs0 = wavread("demo.wav")       # load input signal
	x0 = x0[:,1]
	x0 = resample(x0,Int(Fs)//Int(Fs0)) 
	normalize!(x0,Inf)
	Nt = length(x0)                     # signal length
	x0 .*= 0.8

	C = 0.02                                       # clip threshold
	xc = copy(x0)
	xc[abs.(xc).>=C] = C.*sign.(xc[abs.(xc).>=C])  # clip original signal

	#time win options
	K = 2 # K = 2 hanning
	win = sqrt.(tukey(Nl+1,2/K)[1:Nl])
	overlap = div(Nl,K)

	xd = zeros(x0)                    # allocate output
	x, y = Variable(Nl), Variable(Nl) # optimization variables
	cf = ls( idct(x) - y )            # cost function
	yw = zeros(Nl)

	return Fs, x, x0, xd, xc, y, yw, Nl, Nt, C, win, overlap, cf 
end

function run_demo()
	slv = PANOC(tol = 1e-5, verbose = 0)
	setup = set_up()
	@time solve_problem!(slv, setup..., true)

	Fs, x, x0, xd, xc, y, yw, Nl, Nt, C, overlap, cf = setup
	Irp = find(xc.>=C)  #positive clipping indices
	Irn = find(xc.<=-C) #negative clipping indices 

	SDRy  = 20*log10(norm(x0[[Irp;Irn]])/norm(x0[[Irp;Irn]]-xc[[Irp;Irn]]))
	SDRyd = 20*log10(norm(x0[[Irp;Irn]])/norm(x0[[Irp;Irn]]-xd[[Irp;Irn]]))

	println("SDRy: $(SDRy) SDRyd: $(SDRyd) diff: $(SDRyd-SDRy)")
	return setup
end

function solve_problem!(slv, Fs, x, x0, xd, xc, y, yw, Nl, Nt, C, win, overlap, cf, verb)
	its = Int64[]
	fit_tol = 1e-5
	z = 0                             # frame initial intex
	counter = 0
    winC = C.*win
	#weighted Overlap-Add
	while z+Nl < Nt
		counter += 1
		push!(its,0)

		fill!(x.x,0.) # initialize variables
		fill!(y.x,0.)

		Irp = sort(find(     xc[z+1:z+Nl]  .>=  C) ) #positive clipping indices
		Irn = sort(find(     xc[z+1:z+Nl]  .<= -C) ) #negative clipping indices 
		Im  = sort(find(abs.(xc[z+1:z+Nl]) .<   C)) #non clipped indices
		yw .= xc[z+1:z+Nl].*win

		Mp = GetIndex((Nl,),Irp)
		Mn = GetIndex((Nl,),Irn)
		M  = GetIndex((Nl,),Im )

		if isempty(Irp) && isempty(Irn)
			xd[z+1:z+Nl] .+= yw.*win 
			if verb @printf("%7d / %7d |           no_clip          | \n", z, Nt) end
		else
			N = 0 # number of active components in DCT
			for N = 30:30:30*div(Nl,30)
				cstr = (norm(x,0) <= N, 
					norm(M*y-M*yw) <= sqrt(fit_tol), 
					Mp*y >=  Mp*winC, 
                    Mn*y <= -Mn*winC)
                it, = @minimize cf st cstr with slv
				its[counter] += it
                if norm(idct(~x)- ~y) <= sqrt(fit_tol) break end
			end
			@views xd[z+1:z+Nl] .+= (~y).*win 
			if verb @printf("%7d / %7d | N: %7d /  %7d | \n", z, Nt, N, Nl) end
		end

		z += Nl-overlap                      # update index
	end
	return mean(its)
end

function benchmark(;verb = 0, samples = 1, seconds = 100)

	suite = BenchmarkGroup()

	tol = 1e-4
	solvers = [
               "PANOC",
               "ZeroFPR",
               "PG"
          ]
	slv_opt = ["(verbose = $verb, tol = $tol)", 
		   "(verbose = $verb, tol = $tol)",
		   "(verbose = $verb, tol = $tol)"]

	its = Dict([(sol,0.) for sol in solvers])
	for i in eachindex(solvers)

		setup = set_up()
		solver = eval(parse(solvers[i]*slv_opt[i]))

		suite[solvers[i]] = 
		@benchmarkable(it = solve_problem!(solver, setup..., false), 
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
	println("AudioDecliping (mean of per-frame) its")
	println(its)
	return results
end

function show_results(Fs, x, x0, xd, xc, y, yw, Nl, Nt, C, win, overlap, cf)
	idxs = 2*Nl+1:4*Nl
	t = linspace(0,length(idxs)/Fs,length(idxs)) 
	figure()
	plot(t,x0[idxs])
	plot(t,xd[idxs])
	plot(t,xc[idxs])
end

function save_wav(Fs, x, x0, xd, xc, y, yw, Nl, Nt, C, win, overlap, cf)
    wavwrite(x0 ,"unclipped.wav"; Fs = Fs, nbits = 16, compression=WAVE_FORMAT_PCM)
    wavwrite(xc ,"clipped.wav";   Fs = Fs, nbits = 16, compression=WAVE_FORMAT_PCM)
    wavwrite(0.8*normalize(xd,Inf) ,"declipped.wav"; Fs = Fs, nbits = 16, compression=WAVE_FORMAT_PCM)
end

end
