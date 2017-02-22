
x = OptVar(zeros(5))
A = fft(x)-fft(randn(5))

T1 = 3*norm(x,1)

T2 = ls(A)

cf = ls(A)+3*norm(x,1)
cf = T1+T2

		
println("possible idea for minimize (simplified a lot!)")
println()

main_term = []
idx_main_term = 0
for i in eachindex(cf.Terms)
	if typeof(cf.Terms[i]) <: RegLS.SmoothTerm
		println("search for smooth terms, possibly a least squares with affine operator")
		println()
		show(cf.Terms[i].A)
		main_term = cf.Terms[i].A
		idx_main_term = i
		println()
		println("in case there are other smooth terms I choose the one with the less simple operator")
		println()
		println("otherwise I choose the term with the largest number of variables?")
		println()
	end
end

prox = []
println("check the other terms and obtain their proximal op")
		
println()
for i in eachindex(cf.Terms)
	if i != idx_main_term
		push!(prox, RegLS.get_prox(cf.Terms[i]))
	end
end
show(prox)
		



