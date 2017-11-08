using ProximalOperators
using RegLS
using DSP, WAV

Fs = 16e3                           # sampling frequency
x0, Fs0 = wavread("demo.wav")       # load input signal
x0 = x0[end-Int64(3e5):end]         # cut signal
x0 = resample(x0,Int(Fs)//Int(Fs0)) 
normalize!(x0,Inf)
Nt = length(x0)                     # signal length
x0 .*= 0.8

C = 0.02                                    # clip threshold
xc = copy(x0)
xc[abs.(xc).>=C] = C.*sign.(xc[abs.(xc).>=C])  # clip original signal

Nl = 2^10                     # frames length
#time win options
K = 2 # K = 2 hanning
win = sqrt.(tukey(Nl+1,2/K)[1:Nl])
overlap = div(Nl,K)

xd = zeros(x0)                     # allocate output
x, y = Variable(Nl), Variable(Nl) # optimization variables
cf = ls( idct(x) - y )            # cost function
yw = zeros(Nl)

fit_tol = 1e-5
z = 0                 # frame initial intex

while z+Nl < Nt

	fill!(x.x,0.) # initialize variables
	fill!(y.x,0.)

	Irp = sort(find(    xc[z+1:z+Nl]  .>=  C) ) #positive clipping indices
	Irn = sort(find(    xc[z+1:z+Nl]  .<= -C) ) #negative clipping indices 
	Im  = sort(find(abs.(xc[z+1:z+Nl]) .<   C)) #non clipped indices
	yw .= xc[z+1:z+Nl].*win

	slv = ZeroFPR(tol = 1e-6, verbose = 0)
	N = 0
	for N = 30:30:30*div(Nl,30)
		cstr = (norm(x,0) <= N, 
			norm(y[Im]-yw[Im]) <= fit_tol, 
			y[Irp] in [ C, 1], 
			y[Irn] in [-1,-C])
		slv = @minimize cf st cstr with slv
		if slv.cost <= fit_tol break; end
	end

	@printf("%7d / %7d | N: %7d /  %7d | \n", z, Nt, N, Nl)

	xd[z+1:z+Nl] .+= idct(~x).*win 
	z += Nl-overlap                      # update index
end


using PyPlot
figure()
plot(x0)
plot(xd)
plot(xc)

Irp = find(xc.>=C)  #positive clipping indices
Irn = find(xc.<=-C) #negative clipping indices 

SDRy  = 20*log10(norm(x0[[Irp;Irn]])/norm(x0[[Irp;Irn]]-xc[[Irp;Irn]]))
SDRyd = 20*log10(norm(x0[[Irp;Irn]])/norm(x0[[Irp;Irn]]-xd[[Irp;Irn]]))

println("SDRy: $(SDRy) SDRyd: $(SDRyd) diff: $(SDRyd-SDRy)")


