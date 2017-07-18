__precompile__()

module RegLS

using AbstractOperators, ProximalOperators
import ProximalOperators: RealOrComplex, 
			  is_affine,
			  is_cone,
			  is_convex,
			  is_generalized_quadratic,
			  is_prox_accurate,
			  is_quadratic,
			  is_separable,
			  is_set,
			  is_singleton,
			  is_smooth,
			  is_strongly_convex

include("deep.jl")
include("syntax.jl")
include("solvers.jl")
include("problem.jl")

end
