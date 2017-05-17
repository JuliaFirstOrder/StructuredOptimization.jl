"""
`split(cf::CompositeFunction) -> (smooth, proximable, nonsmooth)`

Splits cost function into `SmoothFunction`, proximable and `NonSmoothFunction` terms.
"""
function split(cf::CompositeFunction)

	smooth,     nonsmooth = split_Smooth(cf)
	proximable, nonsmooth = split_Proximable(nonsmooth)

	return smooth, proximable, nonsmooth

end

"""
`split_Smooth(cf::CompositeFunction) -> (smooth, nonsmooth)`

Splits cost function into `SmoothFunction` and `NonSmoothFunction` terms.
"""
split_Smooth(cf::CompositeFunction) =
CompositeFunction(variable(cf),
	      terms(cf)[is_smooth.(terms(cf))],
	     affine(cf)[is_smooth.(terms(cf))]),
CompositeFunction(variable(cf),
	      terms(cf)[!is_smooth.(terms(cf))],
	     affine(cf)[!is_smooth.(terms(cf))])

"""
`split_Proximable(cf::CompositeFunction) -> (proximable, non_proximable)`

Splits cost function into proximable and non proximable terms.
"""
split_Proximable(cf::CompositeFunction) =
CompositeFunction(variable(cf),
	      terms(cf)[is_gram_diagonal.(operator(cf))],
	     affine(cf)[is_gram_diagonal.(operator(cf))]),
CompositeFunction(variable(cf),
	      terms(cf)[!is_gram_diagonal.(operator(cf))],
	     affine(cf)[!is_gram_diagonal.(operator(cf))])

"""
`split_Quadratic(cf::CompositeFunction) -> (quadratic, non_quadratic)`

Splits cost function into `QuadraticFunction` and non `QuadraticFunction` terms.
"""
split_Quadratic(cf::CompositeFunction) =
CompositeFunction(variable(cf),
	      terms(cf)[is_quadratic.(terms(cf))],
	     affine(cf)[is_quadratic.(terms(cf))]),
CompositeFunction(variable(cf),
	      terms(cf)[!is_quadratic.(terms(cf))],
	     affine(cf)[!is_quadratic.(terms(cf))])
