immutable HCAT <: LinearOperator
	A::AbstractArray{LinearOperator}
	mid::AbstractArray
	function HCAT(A::AbstractArray{LinearOperator}, mid::AbstractArray)
		if any(size.(A[2:end], 1) .!= size(A[1], 1))
			throw(DimensionMismatch("operators must have the same codomain dimension!"))
		end
		if size(mid) != size(A[1], 1)
			throw(DimensionMismatch("buffer must have the correct dimension!"))
		end
		new(A, mid)
	end
end

# So we avoid useless nesting
HCAT(A::LinearOperator) = A

function HCAT(A::Vararg{LinearOperator})
	# we should be able to simplify this
	H = Vector{LinearOperator}()
	mid = zeros(size(A[1],1))
	for a in A push!(H, a) end
	return HCAT(H, mid)
end

size(A::HCAT) = size(A.A[1],1), tuple(size.(A.A, 2)...)

function A_mul_B!(y, S::HCAT, b)
	A_mul_B!(y, S.A[1], b[1])
	for i = 2:length(S.A)
		A_mul_B!(S.mid, S.A[i], b[i])
		y .= (+).(y, S.mid)
	end
end

function At_mul_B!(y, S::HCAT, b)
	for i = 1:length(S.A)
		At_mul_B!(y[i], S.A[i], b)
	end
end

# function *{T<:AbstractArray}(A::HCAT, b::AbstractArray{T, 1})
# 	C = codomainType(A);
# 	y = Array{C}(size(A, 1))
# 	A_mul_B!(y, A, b)
# 	return y
# end

# what is this?
# function .*{T<:AbstractArray}(A::HCAT,b::Array{T,1})
# 	y = Array{AbstractArray,1}(length(A.A))
# 	for i = 1:length(A.A)
# 		y[i] = A.A[i]*b[i]
# 	end
# 	return y
# end

# this is nice, but I don't know if it's needed and correct
# .*(A::HCAT,B::HCAT) = HCAT(A.A.*B.A, A.mid)

import Base: hcat

# define `hcat` for convenience
hcat(A::Vararg{LinearOperator}) = HCAT(A...)

fun_name(S::HCAT) = "Horizontally concatenated operators"

# function fun_dom{D3<:Real}(A::HCAT{D3})
# 	str = ""
# 	for a in A.A str *= fun_D1(a,1) end
# 	str *= "→  ℝ^$(size(A,2))"
# end
#
# function fun_dom{D3<:Complex}(A::HCAT{D3})
# 	str = ""
# 	for a in A.A str *= fun_D1(a,1) end
# 	str *= "→  ℂ^$(size(A,2))"
# end

# fun_D1{D1<:Real, D2}(A::LinearOperator{D1,D2},dim::Int64)    =  " ℝ^$(size(A,dim)) "
# fun_D1{D1<:Complex, D2}(A::LinearOperator{D1,D2},dim::Int64) =  " ℂ^$(size(A,dim)) "

# import Base: copy, sort
#
# copy(A::HCAT) = HCAT(copy(A.A), A.mid)
#
# function sort(A::HCAT,p::Array)
# 	H = A.A[p]
# 	return HCAT(H,A.mid)
# end
#
# extract_operator(A::HCAT, idx::Int64) = A.A[idx]
