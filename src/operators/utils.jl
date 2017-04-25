function optArray{T<:AbstractAffineTerm}(A::T) 
	if typeof(variable(A)) <: Variable 
		return variable(A).x 
	else
		if length(variable(A)) == 1 
			return variable(A)[1].x
		else
			return [variable(A)[i].x for i = 1:length(variable(A))]
		end
	end
end
function optArray!{T<:AbstractAffineTerm,B <:AbstractArray}(A::T,x::B)  
	if typeof(variable(A)) <: Variable 
		copy!(variable(A).x, x)  
	else
		length(variable(A)) != 1 ? error("something went wrong!") : 
		copy!(variable(A)[1].x, x)  
	end
end
function optArray!{T<:AbstractAffineTerm,B <:AbstractArray}(A::T,x::Array{B,1}) 
	for i in eachindex(variable(A))
		copy!(variable(A)[i].x, x[i])  
	end
end

#testing utils
function test_FwAdj(A::LinearOperator, x, y, verb::Bool = false)
	verb && (println(); show(A); println())

	verb && println("forward")
	y = A*x          #verify linear operator works
	verb && @time y = A*x

	y2 = 0*copy(y)
	verb && println("forward preallocated")
	A_mul_B!(y2,A,x) #verify in-place linear operator works
	verb && @time A_mul_B!(y2,A,x)
	test1 =  vecnorm(y-y2) #verify equivalence

	verb && (println(); show(A'); println())
	verb && println("adjoint")
	x = A'*y          #verify adjoint operator inside LinearTerm is the same
	verb && @time x = A'*y

	verb && println("adjoint preallocated")
	x2 = 0*copy(x)
	A_mul_B!(x2,A',y) #verify in-place linear operator works
	verb && @time A_mul_B!(x2,A',y)

	test2 = vecnorm(x-x2) #verify equivalence

	return test1, test2

end

function test_Op(L::LinearOperator,x,y)
	d1 = RegLS.deepvecdot(L*x,  y)
	d2 = RegLS.deepvecdot(x, L'*y)
	return norm( d1 - d2 )   #verify operator and its ajoint
end
