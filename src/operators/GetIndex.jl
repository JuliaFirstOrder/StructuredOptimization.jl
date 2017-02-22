import Base: getindex

immutable GetIndex{D1,D2} <: LinearOp{D1,D2}
	x::OptVar
	idx::Tuple 
	isTranspose::Bool
	dim::Tuple
end
size(A::GetIndex) = A.dim

function getindex{D1,D2}(B::LinearOp{D1,D2}, args...) 
	A = GetIndex{D2,D2}(B.x, args, false, get_size(size(B,2),args...)) 
	return NestedLinearOp(A,B)
end
fun_name(A::GetIndex) = "Get Index"

transpose{D1}(A::GetIndex{D1,D1}) = GetIndex{D1,D1}(A.x, A.idx,true,(A.dim[2],A.dim[1])) 

function A_mul_B!(y::AbstractArray,A::GetIndex,b::AbstractArray) 
	if A.isTranspose
		y.*= 0
		setindex!(y,b,A.idx...) 
	else
		copy!(y,view(b,A.idx...))  
	end
end

function get_size(dim,args...) 
	dim2 = [dim...] 
	for i = 1:length(args)
		if args[i] != Colon() dim2[i] = length(args[i]) end 
	end
	return (dim,tuple(dim2...))
end










