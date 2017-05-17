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

# Constructors

DiagOp{T}(domainType::Type, d::T) =
DiagOp{typeof(d),ndims(d)}(domainType, size(d), d)
DiagOp(d::AbstractArray) = DiagOp{typeof(d),ndims(d)}(eltype(d),size(d), d)
DiagOp{T}(x::T, d::T) = DiagOp{T,ndims(d)}(eltype(x), size(x), d)

# Mappings

function A_mul_B!{T,N}(y::T,L::DiagOp{T,N},b::T)
	y .= (*).(L.d,b)
end

function Ac_mul_B!{T,N}(y::T,L::DiagOp{T,N},b::T)
	y .= (*).(conj.(L.d), b)
end

# Transformations (we'll see about this)
# inv(L::DiagOp) = DiagOp(L.domainType, L.dim_in, (L.d).^(-1))

# Properties

size(L::DiagOp) = (L.dim_in,L.dim_in)

fun_name(L::DiagOp)  = "Diagonal Operator"

is_diagonal(L::DiagOp) = true

# TODO: probably the following allows for too-close-to-singular matrices
is_invertible(L::DiagOp) = all( L.d .!= 0.  )
is_full_row_rank(L::DiagOp) = is_invertible(L)
is_full_column_rank(L::DiagOp) = is_invertible(L)
