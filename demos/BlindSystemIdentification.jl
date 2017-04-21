using RegLS

srand(123)
Fs = 8e3
Nx = 1000    #length of input singal
Nh = 300     #length of RIRs
Nzp = Nx-Nh
SNR = 30

Km = 3 #num of mics

h_true = randn(Nh,Km) # true RIRs

y = zeros(Nx,Km) # preallocate receiver signal
x = randn(Nx)    # input signal 

for i = 1:Km 
	y[:,i] .= filt([h_true[:,i];zeros(Nzp)],[1.],x) # receiver signals
	y[:,i] .+= 10^(-SNR/10)*var(y[:,i])*randn(Nx)   # add noise
end

h = Variable(randn(Nh,Km))
cf = emptycostfun() #empty cost function
	
#create cost function
for i = 1:Km-1,j = i+1:Km
	cf += ls( filt( zeropad(h[:,j],Nzp), y[:,i]) - filt( zeropad(h[:,i],Nzp), y[:,j]) )
end

tol = 1e-8
verb = 0
slv = ZeroFPR(tol = tol, verbose = verb)
#slv = FPG(tol = tol, verbose = verb)

slv = minimize(cf,norm(h) == 1, slv )
println(slv)

# Normalized Projection Misalignment
NPM = 10*log10( vecnorm(h_true-vecdot(~h,h_true)/vecdot(~h,~h)*~h   )^2/vecdot(h_true,h_true) )



