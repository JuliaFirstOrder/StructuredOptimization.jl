
function test_FwAdj(A::LinearOp, x, y)
	println(); show(A); println()

	println("forward")
	y = A*x          #verify linear operator works
	@time y = A*x

	y2 = 0*copy(y)
	println("forward preallocated")
	A_mul_B!(y2,A,x) #verify in-place linear operator works
	@time A_mul_B!(y2,A,x)
	test1 =  vecnorm(y-y2) #verify equivalence

	println("adjoint")
	At = A'
	x = At*y          #verify adjoint operator works
	@time x = At*y

	println("adjoint preallocated")
	x2 = 0*copy(x)
	A_mul_B!(x2,At,y) #verify in-place linear operator works
	@time A_mul_B!(x2,At,y)

	test2 = vecnorm(x-x2) #verify equivalence

	return test1, test2

end

function test_Op(A::LinearOp,x,y)
	return norm( RegLS.deepvecdot(A*x,y) - RegLS.deepvecdot(x,A'*y))   #verify operator and its ajoint
end
