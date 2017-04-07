immutable VCAT <: LinearOperator
	A::AbstractArray{LinearOperator}
	function VCAT(A::AbstractArray{LinearOperator})
		if any(size.(A[2:end], 2) .!= size(A[1], 2))
			throw(DimensionMismatch("operators must have the same domain dimension!"))
		end
		new(A)
	end
end

# So we avoid useless nesting
VCAT(A::LinearOperator) = A

function VCAT(A::Vararg{LinearOperator})
	H = Vector{LinearOperator}()
	for a in A
		push!(H,a)
	end
	return VCAT(H)
end

size(A::VCAT) = tuple(size.(A.A, 1)...), size(A.A[1], 2)

function A_mul_B!(y::AbstractArray, S::VCAT, b::AbstractArray)
	for i = 1:length(S.A)
		A_mul_B!(y[i],S.A[i],b)
	end
end

# function *(A::VCAT,b::AbstractArray)
# 	y = Array{AbstractArray,1}(length(A.A))
# 	for i = 1:length(A.A)
# 		C = codomainType(A.A[i])
# 		y[i] = Array{C}(size(A.A[i],1))
# 	end
# 	A_mul_B!(y,A,b)
# 	return y
# end

# I'm not sure what this is
# .*{D3}(A::VCAT{D3},B::VCAT{D3}) = VCAT{D3}(A.A.*B.A, A.mid)

import Base: vcat

# define `vcat` for convenience
vcat(A::Vararg{LinearOperator}) = VCAT(A...)

fun_name(L::VCAT) = "Vertically concatenated operators"

# #printing stuff
# function fun_dom{D3<:Real}(A::VCAT{D3})
# 	str = " ℝ^$(size(A,1)) → "
# 	for a in A.A str *= fun_D1(a,2) end
# 	return str
# end
#
# function fun_dom{D3<:Complex}(A::VCAT{D3})
# 	str = " ℂ^$(size(A,1)) → "
# 	for a in A.A str *= fun_D1(a,2) end
# 	return str
# end
