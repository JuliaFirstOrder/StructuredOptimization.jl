# StructuredOptimization.jl

`StructuredOptimization.jl` is a high-level modeling language
that utilizes a syntax that is very close to
the mathematical formulation of an optimization problem.

This user-friendly interface
acts as a parser to utilize
three different packages:

* [`ProximalOperators.jl`](https://github.com/kul-forbes/ProximalOperators.jl)

* [`AbstractOperators.jl`](https://github.com/kul-forbes/ProximalOperators.jl)

* [`ProximalAlgorithms.jl`](https://github.com/kul-forbes/ProximalAlgorithms.jl)

`StructuredOptimization.jl` can handle large-scale convex and nonconvex problems with nonsmooth cost functions. 

It supports complex variables as well.

## Installation

From the Julia command line hit `Pkg.clone("https://github.com/nantonel/StructuredOptimization.jl.git")`.
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

See the [documentation]() for more details about the type of problems `StructuredOptimization.jl` can handle and the [demos]() to check out some examples.
