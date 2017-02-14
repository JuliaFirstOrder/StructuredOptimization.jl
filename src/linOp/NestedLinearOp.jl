immutable NestedLinearOp{D1,D2} <: LinearOp{D1,D2}
	x::OptVar{D1}
	y::OptVar{D2}
	A::LinearOp
	B::LinearOp
	status::Array{Int64,1}  
end

function NestedLinearOp(f::Function,B::LinearOp, args...) 
	A = f(B.y, args...)
	return NestedLinearOp(A,B)
end

# constructor with array
function NestedLinearOp(Op::Array) 
	N = NestedLinearOp(Op[end-1],Op[end])
	for i in length(Op)-2:-1:1
		N = NestedLinearOp(Op[i],N)
	end
	return N
end

NestedLinearOp(A::LinearOp,B::LinearOp) = NestedLinearOp(B.x,A.y,A,B,[0])	

function NestedLinearOp(A::LinearOp,B::NestedLinearOp)  
	B.status[1] == 0 ? B.status[1] = 1 : B.status[1] = 2 
		
	NestedLinearOp(B.x,A.y,A,B,[3])	
end

function transpose(N::NestedLinearOp)  
	if N.status == [0]
		NestedLinearOp(N.B', N.A')
	else
		Op = get_op(N)
		Op = flipdim(Op.'[:],1)
		return NestedLinearOp(Op) 
	end
end


*(N::NestedLinearOp,b::AbstractArray) = N.A*(N.B*b) 

function A_mul_B!(y::AbstractArray,N::NestedLinearOp,b::AbstractArray) 
	if N.status[1] == 0
		#println("0")
		A_mul_B!(N.B.y.x, N.B, b      )
		A_mul_B!(y,       N.A, N.B.y.x)
	elseif N.status[1] == 1
		#println("1")
		A_mul_B!(N.B.y.x, N.B, b     )
		A_mul_B!(N.A.y.x, N.A, N.B.y.x)
	elseif N.status[1] == 2
		#println("2")
		A_mul_B!(N.B.y.x, N.B, b      )
		A_mul_B!(N.A.y.x, N.A, N.B.y.x)
	elseif N.status[1] == 3
		#println("3")
		A_mul_B!(N.B.y.x, N.B, b      )
		A_mul_B!(y      , N.A, N.B.y.x)
	end
end

fun_name(N::NestedLinearOp) = N.status[1] == 0 ? string(typeof(N.A))*" * "*string(typeof(N.B)) : "Nested Linear Operator"


#count the number of operators
function number_op(N::NestedLinearOp)
	counter = 0
	if N.status != [0]
		AA = N
		while AA.status != [1]
			counter += 1
			AA = AA.B
		end
	end
	counter += 2
	return counter
end

#get a vector containing all operators
function get_op(N::NestedLinearOp)
	NOp = number_op(N)
	Op = Vector(NOp)
	AA = N
	for i = 1:NOp-2
		Op[i] = AA.A
		AA = AA.B
	end
	Op[end-1] = AA.A
	Op[end]   = AA.B 
	return Op
end
