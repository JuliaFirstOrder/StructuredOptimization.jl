# RegLS.jl

Convex and nonconvex regularized least squares in Julia.

## Installation

From the Julia command line hit `Pkg.clone("https://github.com/nantonel/RegLS.jl.git")`.
Once the package is installed you can update it along with the others issuing
`Pkg.update()` in the command line.

## Usage

With RegLS.jl you can solve problems of the form

```
minimize (1/2)*||L(x) - b||^2 + g(x)
```

Here `L` is a linear operator, `b` is an `Array` of data, and `g` is a regularization function (see next section).
You can use any `AbstractMatrix` object to describe `L`, or any matrix-like object
implementing the matrix-vector product and transpose operations
(see for example [LinearOperators.jl](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl)).
Alternatively, you can provide the direct and adjoint mappings in the form of `Function` objects.

```julia
x, info = solve(L, b, g) # L is a matrix-like object
x, info = solve(Op, OpAdj, b, g, x0) # Op and OpAdj are of type Function
```

The dimensions of `b` must match the ones of `L` or `Op` and `OpAdj`.
Argument `x0` (the initial iterate for the algorithm) is compulsory when
mappings `Op` and `OpAdj` are provided.

## Regularizers

The regularization functions included in `RegLS` are listed here.

Function        | Description                                          | Properties
----------------|------------------------------------------------------|----------------
`ElasticNet`    | Elastic-net regularization                           | convex
`IndAffine`     | Indicator of an affine subspace                      | convex
`IndBallInf`    | Indicator of an infinity-norm ball                   | convex
`IndBallL0`     | Indicator of an L0 pseudo-norm ball                  | nonconvex
`IndBallL2`     | Indicator of an Euclidean ball                       | convex
`IndBallRank`   | Indicator of the set of matrices with given rank     | nonconvex
`IndBox`        | Indicator of a box                                   | convex
`IndHalfspace`  | Indicator of a halfspace                             | convex
`IndNonnegative`| Indicator of the nonnegative orthant                 | convex
`IndSimplex`    | Indicator of the probability simplex                 | convex
`IndSOC`        | Indicator of the second-order cone                   | convex
`NormL0`        | L0 pseudo-norm                                       | nonconvex
`NormL1`        | L1 norm                                              | convex
`NormL2`        | Euclidean norm                                       | convex
`NormL21`       | Sum-of-L2 norms                                      | convex
`SqrNormL2`     | Squared Euclidean norm                               | convex

Each function can be customized with parameters: you can access the documentation of each function from the command line of Julia directly (try typing in `?normL1`).
Once a function has been created, you can at any time inspect it by simply printing it out:

```
julia> g = normL1(3.5);
julia> g
description : weighted L1 norm
type        : C^n → R
expression  : x ↦ λ||x||_1
parameters  : λ = 3.5
```

## Example: sparse signal reconstruction

Consider the problem of recovering a sparse signal, observed through a measurement
matrix `A` which is orthogonal. A random instance of such problem is generated as
follows, where measurements are affected by Gaussian noise:

```julia
m, n, k = 1024, 4096, 160 # problem parameters
B = randn(m, n)
Q, R = qr(B')
A = Q' # measurement matrix
x_orig = sign(randn(n))
J = randperm(n)
x_orig[J[k+1:end]] = 0 # generate sparse -1/+1 signal
sigma = 1e-2
y = A*x_orig + sigma*randn(m) # add noise to measurement
```

One way to approximately reconstruct `x_orig` is to solve an L1-regularized
least squares problem, as in the following snippet (parameters here are taken
	from [this paper](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4407767)):

```julia
using RegLS
lambda_max = norm(A'*y, Inf)
lambda = 0.01*lambda_max # regularization parameter
x_L1, info = solve(A, y, normL1(lambda), zeros(n))
```

Alternatively, one can use the `indBallL0` regularizer to look for the best
approximation with a given number of nonzero coefficients:

```julia
x_L0c, info = solve(A, y, indBallL0(200), zeros(n))
```

## References

The algorithms implemented in RegLS.jl are described in the following papers.

1. L. Stella, A. Themelis, P. Patrinos, “Forward-backward quasi-Newton methods for nonsmooth optimization problems,” [arXiv:1604.08096](http://arxiv.org/abs/1604.08096) (2016).

2. A. Themelis, L. Stella, P. Patrinos, “Forward-backward envelope for the sum of two nonconvex functions: Further properties and nonmonotone line-search algorithms,” [arXiv:1606.06256](http://arxiv.org/abs/1606.06256) (2016).

## Credits

RegLS.jl is developed by [Lorenzo Stella](https://lostella.github.io) and [Niccolò Antonello](http://homes.esat.kuleuven.be/~nantonel/) at [KU Leuven, ESAT/Stadius](https://www.esat.kuleuven.be/stadius/).
