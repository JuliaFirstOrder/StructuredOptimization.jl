function optArray{T<:AffineOperator}(A::T) 
	if typeof(variable(A)) <: OptVar 
		return variable(A).x 
	else
		if length(variable(A)) == 1 
			return variable(A)[1].x
		else
			return [variable(A)[i].x for i = 1:length(variable(A))]
		end
	end
end
function optArray!{T<:AffineOperator,B <:AbstractArray}(A::T,x::B)  
	if typeof(variable(A)) <: OptVar 
		copy!(variable(A).x, x)  
	else
		length(variable(A)) != 1 ? error("something went wrong!") : 
		copy!(variable(A)[1].x, x)  
	end
end
function optArray!{T<:AffineOperator,B <:AbstractArray}(A::T,x::Array{B,1}) 
	for i in eachindex(variable(A))
		copy!(variable(A)[i].x, x[i])  
	end
end

#testing utils
function test_FwAdj(A::LinearOperator, x, y)
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

function test_Op(A::LinearOperator,x,y)
	return norm( RegLS.deepvecdot(A*x,y) - RegLS.deepvecdot(x,A'*y))   #verify operator and its ajoint
end
