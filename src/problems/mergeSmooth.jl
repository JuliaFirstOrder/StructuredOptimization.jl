function mergeSmooth{T<:QuadraticTerm}(x::Array{OptVar,1}, quadratic::Array{T,1})
	if length(quadratic) == 1
		sort!(quadratic[1].A)
		return quadratic[1]
	end
end
