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
  # include("test_proxstuff.jl") # CHECKED
end

@testset "Syntax" begin
  # include("test_variables.jl") # CHECKED
  # include("test_expressions.jl") # CHECKED
  # include("test_AbstractOp_binding.jl") # CHECKED
  # include("test_terms.jl") # CHECKED
end

@testset "Problem construction" begin
  # include("test_problem.jl") # CHECKED
  # include("test_build_minimize.jl") # CHECKED
end

@testset "End-to-end tests" begin
  # include("test_usage_small.jl") # CHECKED
  include("test_usage.jl") # STUFF TO ADJUST
end

end
