#test single block of variable

x = OptVar(zeros(5))
#minimize(ls(fft(x)-fft(randn(5)))+1e-2*norm(x,1) )
#
#x = OptVar(zeros(5))
#minimize(ls(fft(x)-fft(randn(5))), norm(x,1) <= 1/1e-2 )
#
###test 2 blocks of variables
#x, y = OptVar(5), OptVar(Complex{Float64},5)
#minimize(ls(fft(x)+y-fft(randn(5)))+1e-2*norm(x,1) )
#
#x, y = OptVar(5), OptVar(5)
#b = randn(5) 
#X, = minimize(ls(dct(x)+eye(y))+1e-2*norm(x,1), norm(y-b,1) <= 1e-2  )
#
#@test (norm(X[2]-b,1)-1e-2)<1e-8

x, y = OptVar(5), OptVar(5)
b = randn(5) 
X, = minimize(ls(dct(x)+eye(y)), norm(x,2)<=1e2, norm(diagop(y,2)-b,1) <= 1e-2  )
X, = minimize(ls(dct(x)+eye(y)), [norm(x,2)<=1e2, norm(diagop(y,2)-b,1) <= 1e-2] )

@test (norm(2.*X[2]-b,1)-1e-2)<1e-8
@test  norm(X[1],2)<= 1e2

#x = OptVar(zeros(5))
#cf = 1e-8*ls(x)+ls(fft(x)-fft(randn(5)))+1e-2*norm(x,1)
#s,n = RegLS.split(cf)
#show(length(s))



		
#println("possible idea for minimize (simplified a lot!)")
#println()
#
#main_term = []
#idx_main_term = 0
#for i in eachindex(cf.Terms)
#	if typeof(cf.Terms[i]) <: RegLS.SmoothTerm
#		println("search for smooth terms, possibly a least squares with affine operator")
#		println()
#		show(cf.Terms[i].A)
#		main_term = cf.Terms[i].A
#		idx_main_term = i
#		println()
#		println("in case there are other smooth terms I choose the one with the less simple operator")
#		println()
#		println("otherwise I choose the term with the largest number of variables?")
#		println()
#	end
#end
#
#prox = []
#println("check the other terms and obtain their proximal op")
#		
#println()
#for i in eachindex(cf.Terms)
#	if i != idx_main_term
#		push!(prox, RegLS.get_prox(cf.Terms[i]))
#	end
#end
#show(prox)
		



