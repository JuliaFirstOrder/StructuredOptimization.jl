export DiagOp

immutable DiagOp <: DiagonalOperator
	domainType::Type
	dim_in::Tuple
	d::AbstractArray
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
DiagOp(domainType::Type, d::AbstractArray) = DiagOp(domainType, size(d), d)
DiagOp(d::AbstractArray) = DiagOp(eltype(d),size(d), d)
DiagOp(x::AbstractArray, d::AbstractArray) = DiagOp(eltype(x), size(x), d)

# Operators
function A_mul_B!(y::AbstractArray,L::DiagOp,b::AbstractArray)
	y .= (*).(L.d,b)
end

function Ac_mul_B!(y::AbstractArray,L::DiagOp,b::AbstractArray)
	y .= (*).(conj.(L.d),b)
end

# Transformations
inv(L::DiagOp) = DiagOp(L.domainType, L.dim_in, (L.d).^(-1))

# Properties
fun_name(L::DiagOp)  = "Diagonal Operator"

isInvertible(A::DiagOp) = true
