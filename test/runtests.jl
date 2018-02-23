using StructuredOptimization
using AbstractOperators
using AbstractOperators.BlockArrays
using ProximalOperators
using Base.Test
using Base.Profile

srand(0)

@testset "StructuredOptimization" begin

@testset "Calculus" begin
  include("test_proxstuff.jl")
end

@testset "Syntax" begin
  include("test_variables.jl")
  include("test_expressions.jl")
  include("test_AbstractOp_binding.jl")
  include("test_terms.jl")
end

@testset "Problem construction" begin
  include("test_problem.jl")
end

@testset "Integration tests" begin
  include("test_usage_small.jl")
  include("test_usage.jl")
end

end
