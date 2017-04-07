"""
`split(cf::CostFunction) -> (smooth, proximable, nonsmooth)`
 
Splits cost function into `SmoothFunction`, proximable and `NonSmoothFunction` terms.
"""
function split(cf::CostFunction)

	smooth,     nonsmooth = split_Smooth(cf)
	proximable, nonsmooth = split_Proximable(nonsmooth)

	return smooth, proximable, nonsmooth 

end

"""
`split_Smooth(cf::CostFunction) -> (smooth, nonsmooth)`
 
Splits cost function into `SmoothFunction` and `NonSmoothFunction` terms.
"""
split_Smooth(cf::CostFunction) = 
CostFunction(variable(cf), 
	      terms(cf)[isSmooth.(terms(cf))], 
	     affine(cf)[isSmooth.(terms(cf))]),
CostFunction(variable(cf), 
	      terms(cf)[!isSmooth.(terms(cf))], 
	     affine(cf)[!isSmooth.(terms(cf))])

"""
`split_Proximable(cf::CostFunction) -> (proximable, non_proximable)`
 
Splits cost function into proximable and non proximable terms.
"""
split_Proximable(cf::CostFunction) = 
CostFunction(variable(cf), 
	      terms(cf)[isGramDiagonal.(operator(cf))], 
	     affine(cf)[isGramDiagonal.(operator(cf))]),
CostFunction(variable(cf), 
	      terms(cf)[!isGramDiagonal.(operator(cf))], 
	     affine(cf)[!isGramDiagonal.(operator(cf))])

"""
`split_Quadratic(cf::CostFunction) -> (quadratic, non_quadratic)`
 
Splits cost function into `QuadraticFunction` and non `QuadraticFunction` terms.
"""
split_Quadratic(cf::CostFunction) = 
CostFunction(variable(cf), 
	      terms(cf)[isQuadratic.(terms(cf))], 
	     affine(cf)[isQuadratic.(terms(cf))]),
CostFunction(variable(cf), 
	      terms(cf)[!isQuadratic.(terms(cf))], 
	     affine(cf)[!isQuadratic.(terms(cf))])


