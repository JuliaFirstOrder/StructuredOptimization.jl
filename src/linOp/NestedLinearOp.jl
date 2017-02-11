immutable NestedLinearOp <: LinearOp
	x::OptVar
	A::LinearOp
	B::LinearOp
	dim::Tuple
	dom::Tuple
	NestedLinearOp(A::LinearOp,B::LinearOp) = new(B.x,A,B,(B.dim[1],A.dim[2]),(B.dom[1],A.dom[2]))
	function NestedLinearOp(f::Function,B::LinearOp, args...) 
		A = f( OptVar(  reshape(B.x.x,B.dim[2]) ), args... )
		return NestedLinearOp(A,B)
	end
end

transpose(N::NestedLinearOp) = NestedLinearOp(N.B',N.A') 
*(N::NestedLinearOp,b::AbstractArray) = N.A*(N.B*b) 

#TODO improve this?
A_mul_B!(y::AbstractArray,N::NestedLinearOp,b::AbstractArray) = A_mul_B!(y,N.A,N.B*b) 

#fun_name(A::NestedLinearOp) = "n/a"
