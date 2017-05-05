export DiagOp

immutable DiagOp{T,N} <: LinearOperator
	domainType::Type
	dim_in::NTuple{N,Int}
	d::AbstractArray{T,N}
	function DiagOp(domainType,dim_in,d)
		if domainType <: Real && eltype(d) <: Complex
			throw(DomainError())
		end
		if size(d) != dim_in
			throw(DimensionMismatch())
		end
		new(domainType,dim_in,d)
	end
end

size(L::DiagOp) = (L.dim_in,L.dim_in)

# Constructors
DiagOp{T,N}(domainType::Type, d::AbstractArray{T,N}) = 
DiagOp{domainType,N}(domainType, size(d), d)
DiagOp(d::AbstractArray) = DiagOp{eltype(d),ndims(d)}(eltype(d),size(d), d)
DiagOp{T,N}(x::AbstractArray{T,N}, d::AbstractArray{T,N}) = DiagOp{T,N}(eltype(x), size(x), d)

# Operators
function A_mul_B!{T,N}(y::AbstractArray{T,N},L::DiagOp{T,N},b::AbstractArray{T,N})
	y .= (*).(L.d,b)
end

function Ac_mul_B!{T,N}(y::AbstractArray{T,N},L::DiagOp{T,N},b::AbstractArray{T,N})
	y .= (*).(conj.(L.d), b) 
end

# Transformations
inv(L::DiagOp) = DiagOp(L.domainType, L.dim_in, (L.d).^(-1))

# Properties
fun_name(L::DiagOp)  = "Diagonal Operator"

isInvertible(L::DiagOp)    = true
isDiagonal(L::DiagOp)      = true
isGramDiagonal(L::DiagOp)  = true
isFullRowRank(L::DiagOp)   = all( L.d .!= 0.  )
isColumnRowRank(L::DiagOp) = all( L.d .!= 0.  )
