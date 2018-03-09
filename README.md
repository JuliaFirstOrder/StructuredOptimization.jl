# StructuredOptimization.jl

[![Build Status](https://travis-ci.org/kul-forbes/StructuredOptimization.jl.svg?branch=master)](https://travis-ci.org/kul-forbes/StructuredOptimization.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/lfrmkg2s1awyxtk8/branch/master?svg=true)](https://ci.appveyor.com/project/nantonel/abstractoperators-jl/branch/master)
[![codecov](https://codecov.io/gh/kul-forbes/StructuredOptimization.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kul-forbes/StructuredOptimization.jl)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://kul-forbes.github.io/StructuredOptimization.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://kul-forbes.github.io/StructuredOptimization.jl/latest)

`StructuredOptimization.jl` is a high-level modeling language
that utilizes a syntax that is very close to
the mathematical formulation of an optimization problem.

This user-friendly interface
acts as a parser to utilize
three different packages:

* [`ProximalOperators.jl`](https://github.com/kul-forbes/ProximalOperators.jl)

* [`AbstractOperators.jl`](https://github.com/kul-forbes/AbstractOperators.jl)

* [`ProximalAlgorithms.jl`](https://github.com/kul-forbes/ProximalAlgorithms.jl)

`StructuredOptimization.jl` can handle large-scale convex and nonconvex problems with nonsmooth cost functions. 

It supports complex variables as well.

## Installation

From the Julia command line hit `Pkg.clone("https://github.com/kul-forbes/StructuredOptimization.jl.git")`.
Once the package is installed you can update it along with the others issuing `Pkg.update()` in the command line.

## Usage

A *least absolute shrinkage and selection operator* (LASSO) can be solved with only few lines of code:

```julia
julia> using StructuredOptimization

julia> n, m = 100, 10;                # define problem size 

julia> A, y = randn(m,n), randn(m);   # random problem data

julia> x = Variable(n);               # initialize optimization variable

julia> λ = 1e-2*norm(A'*y,Inf);       # define λ    

julia> @minimize ls( A*x - y ) + λ*norm(x, 1); # minimize problem

```

See the [documentation]() for more details about the type of problems `StructuredOptimization.jl` can handle and the [demos](https://github.com/kul-forbes/StructuredOptimization.jl/tree/master/demos) to check out some examples.
