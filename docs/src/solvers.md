# Solvers

## Minimizing a function

```@docs
@minimize
```

!!! note "Problem warm-starting"

    By default *warm-starting* is always enabled.
    For example, if two problems that involve the same variables are solved consecutively,
    the second one will be automatically warm-started by the solution of the first one.
    That is because the variables are always linked to their respective data vectors.
    If one wants to avoid this, the optimization variables needs to be manually re-initialized
    before solving the second problem e.g. to a vector of zeros: `~x .= 0.0`.


## Specifying solver and options

You can pick the algorithm to use as `Solver` object from the
[`ProximalAlgorithms.jl`](https://github.com/kul-forbes/ProximalAlgorithms.jl)
package. Currently, the following algorithms are supported.

```@docs
ForwardBackward
ZeroFPR
PANOC
```


## Parse and solve

The macro [`@minimize`](@ref) automatically parse and solve the problem.
An alternative syntax is given by the function [`problem`](@ref) and [`solve`](@ref).

```@docs
problem
solve
```

Once again, the `Solver` objects is to be picked from
[`ProximalAlgorithms.jl`](https://github.com/kul-forbes/ProximalAlgorithms.jl)).

## References

[[1]](http://www.mit.edu/~dimitrib/PTseng/papers/apgm.pdf) Tseng, *On Accelerated Proximal Gradient Methods for Convex-Concave Optimization* (2008).

[[2]](http://epubs.siam.org/doi/abs/10.1137/080716542) Beck, Teboulle, *A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems*, SIAM Journal on Imaging Sciences, vol. 2, no. 1, pp. 183-202 (2009).

[[3]](https://arxiv.org/abs/1606.06256) Themelis, Stella, Patrinos, *Forward-backward envelope for the sum of two nonconvex functions: Further properties and nonmonotone line-search algorithms*, arXiv:1606.06256 (2016).

[[4]](https://doi.org/10.1109/CDC.2017.8263933) Stella, Themelis, Sopasakis, Patrinos, *A simple and efficient algorithm for nonlinear model predictive control*, 56th IEEE Conference on Decision and Control (2017).
