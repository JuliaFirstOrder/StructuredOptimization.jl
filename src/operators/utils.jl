function optArray{T<:AffineOperator}(A::T) 
	if typeof(A.x) <: OptVar 
		return copy(A.x.x) 
	else
		if length(A.x) == 1 
			return copy(A.x[1].x)
		else
			return [copy(A.x[i].x) for i = 1:length(A.x)]
		end
	end
end
function optArray!{T<:AffineOperator,B <:AbstractArray}(A::T,x::B)  
	if typeof(A.x) <: OptVar 
		copy!(A.x.x, x)  
	else
		length(A.x) != 1 ? error("something went wrong!") : 
		copy!(A.x[1].x, x)  
	end
end
function optArray!{T<:AffineOperator,B <:AbstractArray}(A::T,x::Array{B,1}) 
	for i in eachindex(A.x)
		copy!(A.x[i].x, x[i])  
	end
end

#testing utils
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
