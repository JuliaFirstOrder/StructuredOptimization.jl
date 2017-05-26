export Zeros

immutable Zeros{C, D, O, I} <: LinearOperator
	dim_out::NTuple{O, Integer}
	dim_in::NTuple{I, Integer}
end

# Constructors

Zeros{C, D}(out::Tuple, in::Tuple) where {C, D} =
	Zeros{C, D, length(out), length(in)}(out, in)

# Properties

domainType{C, D, O, I}(L::Zeros{C, D, O, I}) = D
codomainType{C, D, O, I}(L::Zeros{C, D, O, I}) = C

# Mappings

function A_mul_B!{C, D, O, I}(y::AbstractArray{C, O}, A::Zeros{C, D, O, I}, b::AbstractArray{D, I})
	y .= zero(C)
end

function Ac_mul_B!{C, D, O, I}(y::AbstractArray{D, I}, A::Zeros{C, D, O, I}, b::AbstractArray{C, O})
	y .= zero(D)
end

# Properties

size(L::Zeros) = (L.dim_out, L.dim_in)

fun_name(A::Zeros)  = "Zero operator"

is_null(L::Zeros) = true
is_diagonal(L::Zeros) = true
is_gram_diagonal(L::Zeros) = true
