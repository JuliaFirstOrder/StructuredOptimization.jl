#test linear operators

stuff = [
				 Dict( 
					"Operator" => (*,),
					"params" => ((randn(4,4)),),
					"args"   => ( randn(4), randn(4) )
							 ),
				 Dict( 
					"Operator" => (*,),
					"params" => ((randn(4,4)+im*randn(4,4)),),
					"args"   => ( randn(4)+im*randn(4), randn(4)+im*randn(4) )
							 ),
				 Dict( 
					"Operator" => (*,),
					"params" => ((randn(4,6)+im*randn(4,6)),),
					"args"   => ( randn(6)+im*randn(6), randn(4)+im*randn(4) )
							 ),
				 Dict( 
					"Operator" => (fft,),
					"params" => ((), ),
					"args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
							 ),
				 Dict( 
					"Operator" => (fft,),
					"params" => ((),),
					"args"   => ( randn(64,64), randn(64,64)+im*randn(64,64) )
							 ),
				 Dict( 
					"Operator" => (ifft,),
					"params" => ((), ),
					"args"   => ( randn(64,64), randn(64,64)+im*randn(64,64) )
							 ),
				 Dict( 
					"Operator" => (ifft,),
					"params" => ((), ),
					"args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
							 ),
				 Dict( 
					"Operator" => (dct,),
					"params" => ((), ),
					"args"   => ( randn(64,64), randn(64,64) )
							 ),
				 Dict( 
					"Operator" => (idct,),
					"params" => ((), ),
					"args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
							 ),
				 Dict( 
					"Operator" => (reshape,),
					"params" => ((10,10), ),
					"args"   => ( randn(100), randn(10,10) )
							 ),
				 Dict(
					"Operator" => (*, dct),
					"params" => ((randn(32,64)) ,() ),
					"args"   => ( randn(64), randn(32) )
							 ),
				 Dict( 
					"Operator" => (ifft, reshape, dct),
					"params" => ((),(10,10),() ),
					"args"   => ( randn(100), randn(10,10)+im*randn(10,10) )
							 ),
				 Dict(
					"Operator" => (dct, ifft, reshape, dct),
					"params" => ((),(),(10,10),() ),
					"args"   => ( randn(100), randn(10,10)+im*randn(10,10) )
							 ),
				]


for i in eachindex(stuff)

	x,y = deepcopy(stuff[i]["args"])
	X,Y = OptVar(x), OptVar(y)

	params = stuff[i]["params"][1]
	Op     = stuff[i]["Operator"][1]
	(Op == *) ? A = Op(params,X) : A = Op(X, params...)
	for j = 2:length(stuff[i]["params"])
		params = stuff[i]["params"][j]
		Op     = stuff[i]["Operator"][j]
		(Op == *) ? A = Op(params,A) : A = Op(A, params...)
	end
	println(); show(A); println()

	y = A*x          #verify linear operator works
	println("forward")
	@time y = A*x

	y2 = 0*copy(y)
	A_mul_B!(y2,A,x) #verify in-place linear operator works
	println("forward preallocated")
	@time A_mul_B!(y2,A,x)
	@test norm(y-y2) < 1e-8 #verify equivalence

	At = A'
	x = At*y          #verify adjoint operator works
	println("adjoint")
	@time x = At*y

	x2 = 0*copy(x)
	A_mul_B!(x2,At,y) #verify in-place linear operator works
	println("adjoint preallocated")
	@time A_mul_B!(x2,At,y)

	@test norm(x-x2) < 1e-8 #verify equivalence

	X,Y = deepcopy(stuff[i]["args"])
	@test norm( real(vecdot(A*X,Y)) - real(vecdot(X,A'*Y))) <1e-8  #verify operator and its ajoint

end


#x,y = randn(64), randn(32)
#X = OptVar(x)
#A = dct(randn(32,64)*X) 
#
#show(A)
#
#y2 = 0*copy(y)
#y1 = A*x
#@time y1 = A*x
#
#A_mul_B!(y2,A,x)
#@time A_mul_B!(y2,A,x)
#@test norm(y1-y2)<1e-8
#
##Ops,mids = RegLS.disassamble(A)
##show(Ops)
##show(mids)
##N = RegLS.NestedLinearOp(Ops,mids)
##show(N)
#
#x1 = A'*y
#@time x1 = A'*y
#
#x2 = 0*copy(x)
#
#Ac_mul_B!(x2,A,y)
#@time Ac_mul_B!(x2,A,y)
#
#@test norm(x1-x2)<1e-8
#
#x,y = randn(64), randn(32)
#@test norm( (real(vecdot(A*x,y))) - real((vecdot(x,A'*y)))) <1e-8  #verify operator and its ajoint
















