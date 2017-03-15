immutable NestedLinearOperator{D1,D2} <: LinearOperator{D1,D2}
	A::LinearOperator
	B::LinearOperator
	mid::AbstractArray       # memory in the middle of the 2 Op
	dim::Tuple
end
size(A::NestedLinearOperator) = A.dim
variable(A::NestedLinearOperator) = variable(A.B)

function NestedLinearOperator{D1,Dm}(f::Function,B::LinearOperator{D1,Dm}, args...) 
	mid = Array{Dm}(size(B,2))
	(f == *) ? A = f(args[1], OptVar(mid)) : A = f(OptVar(mid), args...)
	return NestedLinearOperator(A, B, mid)
end

*(A::LinearOperator, B::LinearOperator) = NestedLinearOperator(A,B) 
*{E<:IdentityOperator}(A::E, B::LinearOperator) = B

.*(A::LinearOperator, B::LinearOperator) = NestedLinearOperator(A,B) 
.*{E<:IdentityOperator}(A::E, B::LinearOperator) = B

NestedLinearOperator{D1,Dm,D2}(A::LinearOperator{Dm,D2}, B::LinearOperator{D1,Dm}, mid::AbstractArray{Dm}) = 
NestedLinearOperator{D1,D2}(A, B, mid, (size(B,1),size(A,2)))

function NestedLinearOperator{D1,Dm,D2}(A::LinearOperator{Dm,D2},B::LinearOperator{D1,Dm})  
	mid = Array{Dm}(size(B,2))
	return NestedLinearOperator{D1,D2}(A, B, mid, (size(B,1),size(A,2)))	
end


# constructor with array
function NestedLinearOperator(Ops::Array,mids::Array) 
	N = NestedLinearOperator(Ops[end-1],Ops[end],mids[end])
	for i in length(Ops)-2:-1:1
		N = NestedLinearOperator(Ops[i],N,mids[i])
	end
	return N
end

function transpose{D1,D2}(N::NestedLinearOperator{D1,D2})  
	if typeof(N.B) <: NestedLinearOperator
		Ops, mids = disassamble(N)
		Ops = flipdim(Ops.'[:],1)
		mids = flipdim(mids[:],1)
		return NestedLinearOperator(Ops,mids) 
	else
		return NestedLinearOperator(N.B', N.A', N.mid)
	end
end



function A_mul_B!{D1,D2}(y::AbstractArray,N::NestedLinearOperator{D1,D2},b::AbstractArray) 
		A_mul_B!(N.mid,   N.B, b    )
		A_mul_B!(y,       N.A, N.mid)
end

fun_name(N::NestedLinearOperator) = ((typeof(N.B) <: NestedLinearOperator) == false ) ? 
fun_name(N.A)*" * "*fun_name(N.B) : "Nested Linear Operator"


#count the number of operators
function countOp(N::NestedLinearOperator)
	counter = 0
	if typeof(N.B) <: NestedLinearOperator
		AA = N
		while typeof(AA.B) <: NestedLinearOperator
			counter += 1
			AA = AA.B
		end
	end
	counter += 2
	return counter
end
#
##get a vector containing all operators
function disassamble(N::NestedLinearOperator)
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
