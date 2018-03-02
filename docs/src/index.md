# StructuredOptimization.jl

`StructuredOptimization.jl` is a high-level modeling language
that utilizes a syntax that is very close to
the mathematical formulation of an optimization problem.

This user-friendly interface
acts as a parser to utilize
three different packages:

* [`ProximalOperators.jl`](https://github.com/kul-forbes/ProximalOperators.jl) provides proximal mappings of functions that are frequently used in signal processing and optimization. 

* [`AbstractOperators.jl`](https://github.com/kul-forbes/ProximalOperators.jl) provides algorithms for the evaluation and combination of forward and (Jacobian) adjoint of linear and nonlinear mappings.

* [`ProximalAlgorithms.jl`](https://github.com/kul-forbes/ProximalAlgorithms.jl) is a library of proximal algorithms (aka splitting algorithms) solvers.

`StructuredOptimization.jl` can handle large-scale convex and nonconvex problems with nonsmooth cost functions: see ? for a set of demos.

# Credits

`StructuredOptimization.jl` is developed by
[Lorenzo Stella](https://lostella.github.io) and
[Niccolò Antonello](https://nantonel.github.io)
at [KU Leuven, ESAT/Stadius](https://www.esat.kuleuven.be/stadius/).

## Citing

