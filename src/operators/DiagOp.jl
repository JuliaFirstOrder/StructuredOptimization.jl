export DiagOp

immutable DiagOp{A <: AbstractArray, T} <: LinearOperator
	d::A
	# function DiagOp{R, D, T}(d::D{T}) where {R <: Real, D <: AbstractArray, T <: Union{R, Complex{R}}}
	# 	if domainType <: Real && eltype(d) <: Complex
	# 		throw(DomainError())
	# 	end
	# 	new(d, domainType)
	# end
end

# Constructors

DiagOp{A <: AbstractArray}(d::A) = DiagOp{A, eltype(d)}(d)
DiagOp{A <: AbstractArray}(d::A, T::Type) = DiagOp{A, T}(d)

# Mappings

function A_mul_B!{A, T}(y::AbstractArray, L::DiagOp{A, T}, b::AbstractArray)
	y .= (*).(L.d, b)
end

function Ac_mul_B!{A, T}(y::AbstractArray, L::DiagOp{A, T}, b::AbstractArray)
	y .= (*).(conj.(L.d), b)
end

# Transformations (we'll see about this)
# inv(L::DiagOp) = DiagOp(L.domainType, L.dim_in, (L.d).^(-1))

# Properties

domainType{A, T}(L::DiagOp{A, T}) = T
codomainType{A, T}(L::DiagOp{A, T}) = T

size(L::DiagOp) = (size(L.d), size(L.d))

fun_name(L::DiagOp) = "Diagonal Operator"

is_diagonal(L::DiagOp) = true

# TODO: probably the following allows for too-close-to-singular matrices
is_invertible(L::DiagOp) = all(L.d .!= 0.)
is_full_row_rank(L::DiagOp) = is_invertible(L)
is_full_column_rank(L::DiagOp) = is_invertible(L)
