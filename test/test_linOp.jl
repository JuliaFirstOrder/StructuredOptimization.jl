#test linear operators

stuff = [
				 Dict( 
					"Operator" => (eye,),
					"params" => ((),),
					"args"   => ( randn(400), randn(400) )
							 ),
				 Dict( 
					"Operator" => (eye,),
					"params" => ((100),),
					"args"   => ( randn(4), randn(100) )
							 ),
				 Dict( 
					"Operator" => (eye,),
					"params" => ((50),),
					"args"   => ( randn(100), randn(50) )
							 ),
				 Dict( 
					"Operator" => (diagop,),
					"params" => ((randn(2,2)+im*randn(2,2),),),
					"args"   => ( randn(2,2)+im*randn(2,2), randn(2,2)+im*randn(2,2) )
							 ),
				 Dict( 
					"Operator" => (diagop,),
					"params" => ((2,),),
					"args"   => ( randn(2,2), randn(2,2) )
							 ),
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
					"args"   => ( randn(64,64), fft(randn(64,64)) )
							 ),
				 Dict( 
					"Operator" => (ifft,),
					"params" => ((), ),
					"args"   => ( randn(64,64), fft(randn(64,64)) )
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
					"args"   => ( randn(100), dct(reshape(ifft(randn(100)),10,10)) )
							 ),
				 Dict(
					"Operator" => (dct, eye, reshape, dct),
					"params" => ((),(100),(10,10),() ),
					"args"   => ( randn(120), randn(10,10) )
							 ),
				 Dict(
					"Operator" => (ifft, eye),
					"params" => ((),(20) ),
					"args"   => ( randn(200)+im*randn(200), randn(20)+im*randn(20) )
							 ),
				 Dict(
					"Operator" => (fft, getindex),
					"params" => ((),([1:5]) ),
					"args"   => ( randn(12)+im*randn(12), randn(5)+im*randn(5) )
							 ),
				 Dict(
					"Operator" => (dct, getindex),
					"params" => ((),([1:2,:,2:5]) ),
					"args"   => ( randn(5,5,5)+im*randn(5,5,5), randn(2,5,4)+im*randn(2,5,4) )
							 ),
				]


for i in eachindex(stuff)

	x,y = deepcopy(stuff[i]["args"])
	X,Y = OptVar(x), OptVar(y)

	params = stuff[i]["params"][1]
	Op     = stuff[i]["Operator"][1]
	if Op == * 
		A = Op(params,X)
	else
		A = Op(X, params...)
	end
	for j = 2:length(stuff[i]["params"])
		params = stuff[i]["params"][j]
		Op     = stuff[i]["Operator"][j]
		if Op == * 
			A = Op(params,A)
		else
			A = Op(A, params...)
		end
	end
	println(); show(A); println()

	println("forward")
	y = A*x          #verify linear operator works
	@time y = A*x

	y2 = 0*copy(y)
	println("forward preallocated")
	A_mul_B!(y2,A,x) #verify in-place linear operator works
	@time A_mul_B!(y2,A,x)
	@test vecnorm(y-y2) < 1e-8 #verify equivalence

	println("adjoint")
	At = A'
	x = At*y          #verify adjoint operator works
	@time x = At*y

	println("adjoint preallocated")
	x2 = 0*copy(x)
	A_mul_B!(x2,At,y) #verify in-place linear operator works
	@time A_mul_B!(x2,At,y)

	@test vecnorm(x-x2) < 1e-8 #verify equivalence

	X,Y = deepcopy(stuff[i]["args"])
	@test vecnorm( vecdot(A*X,Y) - vecdot(X,A'*Y)) <1e-8  #verify operator and its ajoint

end


# try out stuff

x1,x2,x3 = randn(3),randn(3),randn(3)
X1,X2,X3 = OptVar(x1), OptVar(x2), OptVar(x3)
y = randn(3)
x = [x1,x2,x3]

A = dct(X1)+eye(X1)+dct(X2)+eye(X1)+eye(X3)
show(A)
show(length(A.A))

y = A*x          #verify linear operator works
println("forward")
@time y = A*x

y2 = 0*copy(y)
A_mul_B!(y2,A,x) #verify in-place linear operator works
println("forward preallocated")
@time A_mul_B!(y2,A,x)
@test vecnorm(y-y2) < 1e-8 #verify equivalence

At = A'
x = At*y          #verify adjoint operator works
println("adjoint")
@time x = At*y

x2 = 0*copy(x)
A_mul_B!(x2,At,y) #verify in-place linear operator works
println("adjoint preallocated")
@time A_mul_B!(x2,At,y)

@test vecnorm(x-x2) < 1e-8 #verify equivalence

x1,x2,x3 = randn(3),randn(3),randn(3)
X,Y = [x1,x2,x3],randn(3)

@test norm( (RegLS.deepvecdot(A*X,Y)) - (RegLS.deepvecdot(X,At*Y))) <1e-8  #verify operator and its ajoint


