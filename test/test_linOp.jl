##test linear operators
#
#stuff = [
#				 Dict( "Operator" => fft,
#					     "params" => ((), ),
#							 "args"   => ( randn(64,64), randn(64,64)+im*randn(64,64) )
#							 ),
#				 Dict( "Operator" => fft,
#					     "params" => ((), ),
#							 "args"   => ( randn(64)+im*randn(64), randn(64) )
#							 ),
#				 Dict( "Operator" => fft,
#					     "params" => ((), ),
#							 "args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
#							 ),
#
#				 Dict( "Operator" => ifft,
#					     "params" => ((), ),
#							 "args"   => ( randn(64,64), randn(64,64)+im*randn(64,64) )
#							 ),
#				 Dict( "Operator" => ifft,
#					     "params" => ((), ),
#							 "args"   => ( randn(64)+im*randn(64), randn(64) )
#							 ),
#				 Dict( "Operator" => ifft,
#					     "params" => ((), ),
#							 "args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
#							 ),
#				 Dict( "Operator" => dct,
#					     "params" => ((), ),
#							 "args"   => ( randn(64,64), randn(64,64) )
#							 ),
#				 Dict( "Operator" => idct,
#					     "params" => ((), ),
#							 "args"   => ( randn(64)+im*randn(64), randn(64)+im*randn(64) )
#							 ),
#				 Dict( "Operator" => reshape,
#					     "params" => ((10,10), ),
#							 "args"   => ( randn(100), randn(10,10) )
#							 ),
#
#				]
#
#
#for i in eachindex(stuff)
#	for j in eachindex(stuff[i]["params"])
#
#	x,y = deepcopy(stuff[i]["args"])
#	X,Y = OptVar(x), OptVar(y)
#
#	params = stuff[i]["params"][j]
#	A = stuff[i]["Operator"](X, params...)
#	if j == 1 println(); show(A); println() end
#
#	y = A*x          #verify linear operator works
#	println("forward")
#	@time y = A*x
#
#	y2 = 0*copy(y)
#	A_mul_B!(y2,A,x) #verify in-place linear operator works
#	println("forward preallocated")
#	@time A_mul_B!(y2,A,x)
#	@test norm(y-y2) < 1e-8 #verify equivalence
#
#	x = A'*y          #verify adjoint operator works
#	println("adjoint")
#	@time x = A'*y
#
#	x2 = 0*copy(x)
#	Ac_mul_B!(x2,A,y) #verify in-place linear operator works
#	println("adjoint preallocated")
#	@time Ac_mul_B!(x2,A,y)
#
#	@test norm(x-x2) < 1e-8 #verify equivalence
#
#	X,Y = deepcopy(stuff[i]["args"])
#	@test norm(vecdot(A*X,Y)-vecdot(X,A'*Y)) <1e-8  #verify operator and its ajoint
#end
#
#end

x,y = randn(100,100),randn(10000)
X = OptVar(x)
A = idct(reshape(idct(X),10000))
#A = dct(dct(reshape(X,10000)))
show(A)

y2 = copy(y)
y = A*x
@time y = A*x
A_mul_B!(y2,A,x)
@time A_mul_B!(y2,A,x)

@test norm(y-y2) < 1e-8 #verify equivalence

x2 = copy(x)
x = A'*y
@time x = A'*y
Ac_mul_B!(x2,A,y)
@time Ac_mul_B!(x2,A,y)

@test norm(x-x2) < 1e-8 #verify equivalence


x,y = randn(size(x)),randn(size(y))
@test norm(vecdot(A*x,y)-vecdot(x,A'*y)) <1e-8  #verify operator and its ajoint

