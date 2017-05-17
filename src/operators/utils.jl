#testing utils
function test_FwAdj(A::LinearOperator, x, y, verb::Bool = false)
	verb && (println(); show(A); println())

	verb && println("forward")
	y = A*x          #verify linear operator works
	verb && @time y = A*x

	y2 = 0.*deepcopy(y)
	verb && println("forward preallocated")
	A_mul_B!(y2,A,x) #verify in-place linear operator works
	verb && @time A_mul_B!(y2,A,x)
	test1 =  RegLS.deepvecnorm(y.-y2) #verify equivalence

	verb && (println(); show(A'); println())
	verb && println("adjoint")
	At = A'
	x = At*y          #verify adjoint operator inside LinearTerm is the same
	verb && @time x = At*y

	verb && println("adjoint preallocated")
	x2 = 0.*deepcopy(x)
	A_mul_B!(x2,At,y) #verify in-place linear operator works
	verb && @time A_mul_B!(x2,At,y)

	test2 = RegLS.vecnorm(x.-x2) #verify equivalence

	return test1, test2

end

function test_Op(L::LinearOperator,x,y)
	d1 = RegLS.deepvecdot(L*x,  y)
	d2 = RegLS.deepvecdot(x, L'*y)
	return norm( d1 - d2 )   #verify operator and its ajoint
end
