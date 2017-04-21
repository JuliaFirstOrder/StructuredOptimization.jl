using ProximalOperators
using RegLS
using DSP, WAV

Fs = 16e3
xx, = wavread("demo.wav")
xx = xx[end-Int64(3e5):end]
if Fs == 16e3
	xx = resample(xx,4//11) #for 16 kHz
elseif Fs == 44e3

else 
	error("sampling freq only at 44e3 or 16e3!")
end

xx .*= 0.8/maximum(abs(xx)) #original signal
C = 0.02
y_clip = copy(xx)
y_clip[abs(y_clip).>=C] = C.*sign(y_clip[abs(y_clip).>=C])

Nl = 2^10 #frame length
overlap = round(Int64,Nl*0.5)
y_f = arraysplit(y_clip[:],Nl,overlap) #split array into frames
x_f = arraysplit(xx[:],Nl,overlap)     #split array into frames
y_df = [zeros(Nl) for i = 1:length(y_f) ] #declipped frames
wa = sqrt(hamming(Nl))    #analysis window
ws = sqrt(hamming(Nl))    #synthesis window

x, y = Variable(Nl), Variable(Nl)

cf = ls( idct(x) -y )

fit_tol = 1e-5
greedy = true
reg = IndBallL0
#reg = NormL1
slv = ZeroFPR(tol = 1e-6, verbose = 0)

tic()
for f = 1:length(y_f)
	Irp = find(    y_f[f]  .>=  C)  #positive clipping indices
	Irn = find(    y_f[f]  .<= -C) #negative clipping indices 
	Im  = find(abs(y_f[f]) .<   C) #non clipped indices

	if countnz(Irp)+countnz(Irn)!=0 #check the frame is clipped
		yw = y_f[f].*wa           #windowed frame

		~x .= ones(~x)
		~y .= yw

		k = 0
		if greedy == true
			if reg == IndBallL0
				for k=20:20:20*div(Nl,20)
					cstr = [
	     					norm(x,0)<= k, 
						norm(y[Im]-yw[Im])<=fit_tol, 
						y[Irp] in [ yw[Irp], wa[Irp]], 
						y[Irn] in [-wa[Irn], yw[Irn]] 
					       ]
					slv = minimize(cf, cstr, slv)

					if slv.cost <= fit_tol break; end
				end
			elseif reg == NormL1
				for k in logspace(-1,-5,40)
					#TODO ?
				end
				##
			end
		else
			if reg == IndBallL0
				#TODO ?
			elseif reg == NormL1
				#TODO ?
			end
		end
		@printf("frame: %3d / %3d | nnz: %3d \n ", f, length(y_f), k)
		y_df[f] = idct(~x)

	else #do not process frame
		y_df[f] = y_f[f].*wa  
	end
end
toc()

#overlap and add TODO improve this
y_declip = zeros(xx)
Wnorm = zeros(xx)
counter = 0
k = 0
#x_f[end] = ones(Nl) 
OLA = zeros(Nl)
for k = 1:length(y_df)
	OLA = y_df[k].*ws 
	y_declip[counter+1:counter+Nl] += OLA
	Wnorm[counter+1:counter+Nl] += ws.*wa
	counter += Nl-overlap
end
y_declip ./= Wnorm 
y_declip = y_declip[1:end-Nl]
xx = xx[1:end-Nl]
y_clip = y_clip[1:length(y_declip)]


using PyPlot
figure()
plot(xx)
plot(y_declip)
plot(y_clip)

Irp = find(y_clip.>=C)  #positive clipping indices
Irn = find(y_clip.<=-C) #negative clipping indices 

SDRy  = 20*log10(norm(xx[[Irp;Irn]])/norm(xx[[Irp;Irn]]-y_clip[[Irp;Irn]]))
SDRyd = 20*log10(norm(xx[[Irp;Irn]])/norm(xx[[Irp;Irn]]-y_declip[[Irp;Irn]]))

println("SDRy: $(SDRy) SDRyd: $(SDRyd) diff: $(SDRyd-SDRy)")



