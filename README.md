# JuLiASSO

Welcome to JuLiASSO, a Julia package to solve regularized least squares problems of the form

	minimize (1/2)||Ax-b||^2 + g(x)

where `A` is a linear operator and `g` is a regularization term.

## Installation

From the Julia command line hit:

	Pkg.clone("https://github.com/nantonel/JuLiASSO.jl.git")

Once the package is installed you can update it with:

	Pkg.update()

## Usage

For a matrix `A` you can use:

	using JuLiASSO
	JuLiASSO.solve(A, b, g, x0)

For functions `Op` and `OpAdj` computing the direct and adjoint operator respectively you should call instead:

	using JuLiASSO
	JuLiASSO.solve(Op, OpAdj, b, g, x0)

Argument `x0` is the initial approximation to the solution. Both `x0` and `b` must be `Array{}` objects whose dimensions are conformant with those of `A` or `Op` and `OpAdj`.

## Regularization

Argument `g` in the examples above is the regularization term in the problem. The regularizers included in `JuLiASSO` are listed here.

Function        | Description
----------------|------------------------------------
`indBallInf`    | Indicator of an infinity-norm ball
`indBallL0`     | Indicator of an L0 pseudo-norm ball
`indBallRank`   | Indicator of the set of matrices with given rank
`indBox`        | Indicator of a box
`normL0`        | L0 pseudo-norm
`normL1`        | L1 norm
`normL2`        | Euclidean norm
`normL21`       | Sum-of-L2 norms
`normL2sqr`     | Squared Euclidean norm

Each function can be customized with parameters. You can access the full documentation of each of these functions from the command line of Julia directly:

	julia> ?JuLiASSO.normL1
		normL1(λ::Array{Float64})

		Returns the function g(x) = sum(λ_i|x_i|, i = 1,...,n), for a vector of real
		parameters λ_i ⩾ 0.

		normL1(λ::Float64=1.0)

		Returns the function g(x) = λ||x||_1, for a real parameter λ ⩾ 0.
		
## Example - Some nice example #1

## Example - Some nice example #2

## References
