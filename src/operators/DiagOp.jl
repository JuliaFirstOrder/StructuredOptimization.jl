export DiagOp

immutable DiagOp{T<:AbstractArray,N} <: LinearOperator
	domainType::Type
	dim_in::NTuple{N,Int}
	d::T
	function DiagOp{T,N}(domainType,dim_in,d) where {T<:AbstractArray,N}
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
DiagOp{T}(domainType::Type, d::T) = 
DiagOp{typeof(d),ndims(d)}(domainType, size(d), d)
DiagOp(d::AbstractArray) = DiagOp{typeof(d),ndims(d)}(eltype(d),size(d), d)
DiagOp{T}(x::T, d::T) = DiagOp{T,ndims(d)}(eltype(x), size(x), d)

# Operators
function A_mul_B!{T,N}(y::T,L::DiagOp{T,N},b::T)
	y .= (*).(L.d,b)
end

function Ac_mul_B!{T,N}(y::T,L::DiagOp{T,N},b::T)
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
