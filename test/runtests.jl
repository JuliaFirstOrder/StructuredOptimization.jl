using StructuredOptimization
using AbstractOperators
using ProximalOperators
using Base.Test
using Base.Profile

srand(0)

@testset "StructuredOptimization" begin

@testset "Calculus" begin
  include("test_proxstuff.jl")
end

@testset "Syntax" begin
  include("test_affine.jl")
  include("test_terms.jl")
end

# this test are essentially ProximalAlgorithms tests
# TODO move to ProximalAlgorithms?
#@testset "Solvers" begin  
#  include("test_solvers.jl")
#end

@testset "Problem construction" begin
  include("test_problem.jl")
end

@testset "Integration tests" begin
  include("test_usage_small.jl")
  include("test_usage.jl")
end

end
