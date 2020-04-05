# StructuredOptimization.jl

[![Build status](https://github.com/kul-forbes/StructuredOptimization.jl/workflows/CI/badge.svg)](https://github.com/kul-forbes/StructuredOptimization.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/kul-forbes/StructuredOptimization.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kul-forbes/StructuredOptimization.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kul-forbes.github.io/StructuredOptimization.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kul-forbes.github.io/StructuredOptimization.jl/latest)

StructuredOptimization.jl is a high-level modeling language
that utilizes a syntax that is very close to
the mathematical formulation of an optimization problem.

This user-friendly interface
acts as a parser to utilize
three different packages:

* [ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl)

* [AbstractOperators.jl](https://github.com/kul-forbes/AbstractOperators.jl)

* [ProximalAlgorithms.jl](https://github.com/kul-forbes/ProximalAlgorithms.jl)

StructuredOptimization.jl can handle large-scale convex and nonconvex problems with nonsmooth cost functions.

It supports complex variables as well.

## Installation

To install the package, hit `]` from the Julia command line to enter the package manager, then

```julia
pkg> add StructuredOptimization
```

## Usage

A *least absolute shrinkage and selection operator* (LASSO) can be solved with only few lines of code:

```julia
julia> using StructuredOptimization

julia> n, m = 100, 10;                # define problem size

julia> A, y = randn(m,n), randn(m);   # random problem data

julia> x = Variable(n);               # initialize optimization variable

julia> λ = 1e-2*norm(A'*y,Inf);       # define λ    

julia> @minimize ls( A*x - y ) + λ*norm(x, 1); # solve problem

julia> ~x                             # inspect solution
100-element Array{Float64,1}:
  0.0
  0.0
  0.0
  0.440254
  0.0
  0.0
  0.0
[...]
```

See the [documentation](https://kul-forbes.github.io/StructuredOptimization.jl/latest) for more details about the type of problems StructuredOptimization.jl can handle and the [demos](https://kul-forbes.github.io/StructuredOptimization.jl/stable/demos/) to check out some examples.
