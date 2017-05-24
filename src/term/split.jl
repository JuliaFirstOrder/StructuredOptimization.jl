"""
`split(cf::Term) -> (smooth, proximable, nonsmooth)`

Splits cost function into `SmoothFunction`, proximable and `NonSmoothFunction` terms.
"""
function split(cf::Term)

	smooth,     nonsmooth = split_Smooth(cf)
	proximable, nonsmooth = split_Proximable(nonsmooth)

	return smooth, proximable, nonsmooth

end

"""
`split_Smooth(cf::Term) -> (smooth, nonsmooth)`

Splits cost function into `SmoothFunction` and `NonSmoothFunction` terms.
"""
split_Smooth(cf::Term) =
Term(variable(cf),
	      terms(cf)[is_smooth.(terms(cf))],
	     affine(cf)[is_smooth.(terms(cf))]),
Term(variable(cf),
	      terms(cf)[!is_smooth.(terms(cf))],
	     affine(cf)[!is_smooth.(terms(cf))])

"""
`split_Proximable(cf::Term) -> (proximable, non_proximable)`

Splits cost function into proximable and non proximable terms.
"""
split_Proximable(cf::Term) =
Term(variable(cf),
	      terms(cf)[is_gram_diagonal.(operator(cf))],
	     affine(cf)[is_gram_diagonal.(operator(cf))]),
Term(variable(cf),
	      terms(cf)[!is_gram_diagonal.(operator(cf))],
	     affine(cf)[!is_gram_diagonal.(operator(cf))])

"""
`split_Quadratic(cf::Term) -> (quadratic, non_quadratic)`

Splits cost function into `QuadraticFunction` and non `QuadraticFunction` terms.
"""
split_Quadratic(cf::Term) =
Term(variable(cf),
	      terms(cf)[is_quadratic.(terms(cf))],
	     affine(cf)[is_quadratic.(terms(cf))]),
Term(variable(cf),
	      terms(cf)[!is_quadratic.(terms(cf))],
	     affine(cf)[!is_quadratic.(terms(cf))])
