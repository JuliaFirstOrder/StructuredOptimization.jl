__precompile__()

module RegLS

using AbstractOperators, ProximalOperators
import ProximalOperators: RealOrComplex
import AbstractOperators: domainType, 
			  codomainType,
			  is_eye,
			  is_null,
			  is_diagonal,
			  is_gram_diagonal,
			  is_invertible,
			  is_full_row_rank,
			  is_full_column_rank


include("deep.jl")
include("functions.jl")
include("syntax.jl")
include("solvers.jl")
include("problem.jl")

end
