using RegLS
using AbstractOperators
using ProximalOperators
using Base.Test
using Base.Profile

srand(0)

### Passing:

@testset "RegLS" begin

@testset "Tuple operations" begin
  include("test_deep.jl")
end

@testset "Functions" begin
  include("test_functions.jl")
end

@testset "syntax/AffineExpression" begin
  include("test_variables_expressions.jl")
end

@testset "syntax/Terms" begin
  include("test_terms.jl")
end

@testset "Problems construction" begin
  include("test_problem.jl")
end

@testset "Solvers" begin
  include("test_solvers.jl")
end

end

### Yet to be checked:

# include("test_costfunction.jl")
# include("test_matcomp.jl")
# include("test_svm.jl")
# include("test_minimize.jl")
# include("test_solvers_lasso.jl")
# include("test_solvers.jl")
