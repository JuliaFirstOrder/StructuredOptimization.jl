using RegLS
using AbstractOperators
using ProximalOperators
using Base.Test
using Base.Profile

srand(0)

@testset "RegLS" begin

@testset "Tuple operations" begin
  include("test_deep.jl")
end

@testset "TO MOVE ProximalOperators" begin
  include("test_proxstuff.jl")
end

@testset "syntax/AffineExpression" begin
  include("test_affine.jl")
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

@testset "Usage" begin
  include("test_usage_small.jl")
  include("test_usage.jl")
end

end
