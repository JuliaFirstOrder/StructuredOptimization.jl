
"""
`split_smooth(cf::Vararg{Term}) -> (smooth, nonsmooth)`

Splits cost function into `SmoothFunction` and `NonSmoothFunction` terms.
"""
split_smooth(cf::Vararg{Term}) = cf[find(is_smooth(cf))],cf[find((!).(is_smooth(cf)))]
split_smooth{N}(cf::NTuple{N,Term}) = split_smooth(cf...)

"""
`split_gram_diagonal(cf::Vararg{Term}) -> (proximable, non_proximable)`

Splits cost function into terms with gram diagonal operators and non gram diagonal operators.
"""
split_gram_diagonal(cf::Vararg{Term}) = cf[find(is_gram_diagonal(cf))],cf[find((!).(is_gram_diagonal(cf)))]
split_gram_diagonal{N}(cf::NTuple{N,Term}) = split_gram_diagonal(cf...)

#""" TODO
#`split_Quadratic(cf::Vararg{Term}) -> (quadratic, non_quadratic)`
#
#Splits cost function into `QuadraticFunction` and non `QuadraticFunction` terms.
#"""
