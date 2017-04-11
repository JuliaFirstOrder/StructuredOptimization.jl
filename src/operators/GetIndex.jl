import Base: getindex
export GetIndex 

immutable GetIndex <: LinearOperator
	domainType::Type
	dim_in::Tuple
	dim_out::Tuple
	idx::Tuple

	function GetIndex(domainType,dim_in,idx)
		dim_out = get_dim_out(dim_in,idx...)
		new(domainType,dim_in,dim_out,idx)
	end
end

size(L::GetIndex) = (L.dim_out,L.dim_in)

# Constructors
GetIndex(dim_in::Tuple, idx::Tuple) = GetIndex(Float64, dim_in, idx)
GetIndex(x::AbstractArray, idx::Tuple) = GetIndex(eltype(x), size(x), idx)

getindex(A::LinearOperator,idx...) = GetIndex(codomainType(A),size(A,1),idx)*A 

# Operators
function A_mul_B!(y::AbstractArray,L::GetIndex,b::AbstractArray)
	copy!(y,view(b,L.idx...))
end

function Ac_mul_B!(y::AbstractArray,L::GetIndex,b::AbstractArray)
	y .= 0.
	setindex!(y,b,L.idx...)
end

# Properties
fun_name(L::GetIndex) = "Get Index"
isGramDiagonal(L::GetIndex) = true

# Utils
get_idx(L::GetIndex) = L.idx

function get_dim_out(dim,args...)
	if length(args) != 1
		dim2 = [dim...]
		for i = 1:length(args)
			if args[i] != Colon() dim2[i] = length(args[i]) end
		end
		return tuple(dim2...)
	else
		return tuple(length(args[1]))
	end
end

#################################################################################
## FROM HERE ON IT IS USERS' SYNTAX
#################################################################################
#
#function getindex(x::Variable, args...)
#	A = GetIndex(args, false, get_size(size(x),args...))
#	Affine([x], A, A', Nullable{AbstractArray}() )
#end
#
#function getindex(B::AffineOperator, args...)
#	A = GetIndex{domainType(B),codomainType(B)}(args, false,
#					     get_size(size(operator(B),2),args...))
#	N = NestedLinearOperator(A,operator(B))
#	b = Nullable{AbstractArray}()
#	isnull(B.b) ? nothing : b = adjoint(A)*get(B.b)
#	Affine(variable(B),N,N',b)
#end
