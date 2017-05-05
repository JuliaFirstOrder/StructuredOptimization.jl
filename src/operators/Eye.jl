export Eye

immutable Eye{N} <: LinearOperator
	domainType::Type
	dim_in::NTuple{N,Int}
end

size(L::Eye) = (L.dim_in, L.dim_in)

# Constructors
Eye(dim_in::Vararg{Int64})  = Eye(dim_in)
Eye(domainType::Type, dim_in::Vararg{Int64})  = Eye(domainType,dim_in)
Eye(dim_in::Tuple)          = Eye(Float64,dim_in)
Eye{T}(x::AbstractArray{T}) = Eye(T,size(x))

# Operators
A_mul_B!{T}(y::AbstractArray{T}, L::Eye, b::AbstractArray{T}) = y .= b
Ac_mul_B!(y, L::Eye, b) = A_mul_B!(y,L,b) 

# Transformations
transpose(L::Eye) = L
inv(L::Eye ) = L

#Properties
fun_name(L::Eye)  = "Identity Operator"

           isEye(L::Eye) = true 
      isDiagonal(L::Eye) = true
  isGramDiagonal(L::Eye) = true
    isInvertible(L::Eye) = true
   isFullRowRank(L::Eye) = true
isFullColumnRank(L::Eye) = true

#utils
