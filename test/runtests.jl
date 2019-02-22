using StructuredOptimization
using AbstractOperators
using ProximalOperators
using ProximalAlgorithms
using RecursiveArrayTools
using LinearAlgebra, Random
using DSP, FFTW
using Test

Random.seed!(0)

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
  include("test_build_minimize.jl")
end

@testset "Integration tests" begin
  include("test_usage_small.jl")
  include("test_usage.jl")
end

end
