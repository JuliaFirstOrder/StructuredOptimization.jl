# JuLiASSO

Welcome to JuLiASSO, a Julia package to solve regularized least squares problems of the form

	minimize (1/2)||Ax-b||^2 + g(x)

where `A` is a linear operator and `g` is a regularization term.

## Installation

From the Julia command line hit

	Pkg.clone("https://github.com/nantonel/JuLiASSO.jl.git")
	
Once the package is installed you can update it with

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

Argument `g` in the examples above is the regularization term in the problem. The regularizers included in `JuLiASSO` are listed here

	indBallInf       normL0
	indBallL0        normL1
	IndBallRank      normL2
	indBox           normL21

The names of the available functions are meant to be somehow self-explicative. However, you can access the documentation of each of these from the command line of Julia directly:

	julia> ?
	help?> JuLiASSO.normL1
	  No documentation found.
	
	  JuLiASSO.normL1 is a generic Function.

## Example - Some nice example #1

## Example - Some nice example #2

## References