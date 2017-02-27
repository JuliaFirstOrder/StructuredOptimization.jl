
function Primal(smooth::Array{OptTerm,1}, nonsmooth::Array{OptTerm,1})
	#analyze smooth terms
	if length(smooth) == 0 error("no smooth term is present!") end
	f = smooth[1]
	reg = []
	if length(smooth) > 1
		error("merge regularizers not implemented yet")
	end

	x = f.A.x #extract variables
	g = get_prox(x,nonsmooth)
	 
	return f.A, g
	
end
