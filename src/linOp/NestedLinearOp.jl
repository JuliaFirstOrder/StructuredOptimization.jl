immutable NestedLinearOp{D1,D2} <: LinearOp{D1,D2}
	A::LinearOp
	B::LinearOp
	mid::AbstractArray       # memory in the middle of the 2 Op
	dim::Tuple
end
size(A::NestedLinearOp) = A.dim

function NestedLinearOp{D1,Dm}(f::Function,B::LinearOp{D1,Dm}, args...) 
	mid = Array{Dm}(size(B,2))
	(f == *) ? A = f(args[1], OptVar(mid)) : A = f(OptVar(mid), args...)
	return NestedLinearOp(A,B,mid)
end
#TODO: alternatevely construct using A*B*C...?

NestedLinearOp{D1,Dm,D2}(A::LinearOp{Dm,D2},B::LinearOp{D1,Dm},mid::AbstractArray{Dm}) = 
NestedLinearOp{D1,D2}(A,B,mid,(size(B,1),size(A,2)))

function NestedLinearOp{D1,Dm,D2}(A::LinearOp{Dm,D2},B::LinearOp{D1,Dm})  
	mid = Array{Dm}(size(B,2))
	return NestedLinearOp{D1,D2}(A,B,mid,(size(B,1),size(A,2)))	
end


# constructor with array
function NestedLinearOp(Ops::Array,mids::Array) 
	N = NestedLinearOp(Ops[end-1],Ops[end],mids[end])
	for i in length(Ops)-2:-1:1
		N = NestedLinearOp(Ops[i],N,mids[i])
	end
	return N
end

function transpose{D1,D2}(N::NestedLinearOp{D1,D2})  
	if typeof(N.B) <: NestedLinearOp
		Ops, mids = disassamble(N)
		Ops = flipdim(Ops.'[:],1)
		mids = flipdim(mids[:],1)
		return NestedLinearOp(Ops,mids) 
	else
		return NestedLinearOp(N.B', N.A', N.mid)
	end
end



function A_mul_B!{D1,D2}(y::AbstractArray,N::NestedLinearOp{D1,D2},b::AbstractArray) 
		A_mul_B!(N.mid,   N.B, b    )
		A_mul_B!(y,       N.A, N.mid)
end

fun_name(N::NestedLinearOp) = ((typeof(N.B) <: NestedLinearOp) == false ) ? 
fun_name(N.A)*" * "*fun_name(N.B) : "Nested Linear Operator"


#count the number of operators
function countOp(N::NestedLinearOp)
	counter = 0
	if typeof(N.B) <: NestedLinearOp
		AA = N
		while typeof(AA.B) <: NestedLinearOp
			counter += 1
			AA = AA.B
		end
	end
	counter += 2
	return counter
end
#
##get a vector containing all operators
function disassamble(N::NestedLinearOp)
	NOp = countOp(N)
	Ops   = Vector(NOp)
	mids = Vector(NOp-1)
	AA = N
	for i = 1:NOp-2
		Ops[i] = AA.A
		mids[i] = AA.mid
		AA = AA.B
	end
	mids[end] = AA.mid
	Ops[end-1] = AA.A
	Ops[end]   = AA.B 
	return Ops,mids
end
