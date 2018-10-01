# StructuredOptimization.jl

StructuredOptimization.jl is a high-level modeling language
that utilizes a syntax that is very close to
the mathematical formulation of an optimization problem.

This user-friendly interface
acts as a parser to utilize
three different packages:

* [ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl) provides proximal mappings of functions that are frequently used in signal processing and optimization.

* [AbstractOperators.jl](https://github.com/kul-forbes/AbstractOperators.jl) provides algorithms for the evaluation and combination of forward and (Jacobian) adjoint of linear and nonlinear mappings.

* [ProximalAlgorithms.jl](https://github.com/kul-forbes/ProximalAlgorithms.jl) is a library of proximal algorithms (aka splitting algorithms) solvers.

StructuredOptimization.jl can handle large-scale convex and nonconvex problems with nonsmooth cost functions. It supports complex variables as well. See the [Quick tutorial guide](@ref) and the [Demos](@ref).

## Installation

To install the package, hit `]` from the Julia command line to enter the package manager, then

```julia
pkg> add StructuredOptimization
```

## Citing

If you use StructuredOptimization.jl for published work, we encourage you to cite:

* N. Antonello, L. Stella, P. Patrinos, T. van Waterschoot, “Proximal Gradient Algorithms: Applications in Signal Processing,” [arXiv:1803.01621](https://arxiv.org/abs/1803.01621) (2018).

## Credits

StructuredOptimization.jl is developed by
[Lorenzo Stella](https://lostella.github.io) and
[Niccolò Antonello](https://nantonel.github.io)
at [KU Leuven, ESAT/Stadius](https://www.esat.kuleuven.be/stadius/).
