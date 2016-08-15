# RegLS

Convex and nonconvex regularized least squares in Julia.

## Installation

From the Julia command line hit:

	Pkg.clone("https://github.com/nantonel/RegLS.jl.git")

Once the package is installed you can update it with:

	Pkg.update()

## Usage

For a matrix `A` you can use:

	x = RegLS.solve(A, b, g, x0)

For functions `Op` and `OpAdj` computing the direct and adjoint operator respectively you should call instead:

	x = RegLS.solve(Op, OpAdj, b, g, x0)

Argument `x0` is the initial approximation to the solution. Both `x0` and `b` must be `Array{}` objects whose dimensions are conformant with those of `A` or `Op` and `OpAdj`.

## Regularizers

Argument `g` in the examples above is the regularization term in the problem. The regularizers included in `RegLS` are listed here.

Function        | Description                                          | Properties
----------------|------------------------------------------------------|----------------
`indBallInf`    | Indicator of an infinity-norm ball                   | convex
`indBallL0`     | Indicator of an L0 pseudo-norm ball                  | nonconvex
`indBallRank`   | Indicator of the set of matrices with given rank     | nonconvex
`indBox`        | Indicator of a box                                   | convex
`normL0`        | L0 pseudo-norm                                       | nonconvex
`normL1`        | L1 norm                                              | convex
`normL2`        | Euclidean norm                                       | convex
`normL2sqr`     | Squared Euclidean norm                               | convex
`normL21`       | Sum-of-L2 norms                                      | convex

Each function can be customized with parameters. You can access the full documentation of each of these functions from the command line of Julia directly:

	julia> ?RegLS.normL1
		normL1(λ::Array{Float64})

		Returns the function g(x) = sum(λ_i|x_i|, i = 1,...,n), for a vector of real
		parameters λ_i ⩾ 0.

		normL1(λ::Float64=1.0)

		Returns the function g(x) = λ||x||_1, for a real parameter λ ⩾ 0.

## Example - Some nice example #1

## Example - Some nice example #2

## References

The algorithms implemented in RegLS are described in the following papers.

1. L. Stella, A. Themelis, P. Patrinos, “Forward-backward quasi-Newton methods for nonsmooth optimization problems,” [arXiv:1604.08096](http://arxiv.org/abs/1604.08096) (2016).

2. A. Themelis, L. Stella, P. Patrinos, “Forward-backward envelope for the sum of two nonconvex functions: Further properties and nonmonotone line-search algorithms,” [arXiv:1606.06256](http://arxiv.org/abs/1606.06256) (2016).