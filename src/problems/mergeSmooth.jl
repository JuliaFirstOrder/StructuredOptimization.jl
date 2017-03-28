function mergeSmooth(x_sorted::Array{OptVar,1}, cf::CostFunction)
	sA = Vector{AffineOperator}(length(affine(cf)))
	for i in eachindex(affine(cf)) 
		sA[i] = sort_and_expand(x_sorted,affine(cf)[i])
	end
	return CostFunction(x_sorted,cf.f,sA)
end
