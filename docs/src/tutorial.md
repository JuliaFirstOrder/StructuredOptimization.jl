# Quick tutorial guide

## Standard problem formulation

Currently with StructuredOptimization.jl one can solve problems of the form

```math
\underset{ \mathbf{x} }{\text{minimize}} \ f(\mathbf{x}) + g(\mathbf{x}),
```

where $f$ is a smooth function while $g$ is possibly nonsmooth.

## Unconstrained optimization

The *least absolute shrinkage and selection operator* (LASSO) belongs to this class of problems:

```math
\underset{ \mathbf{x} }{\text{minimize}} \ \tfrac{1}{2} \| \mathbf{A} \mathbf{x} - \mathbf{y} \|^2+ \lambda \| \mathbf{x} \|_1.
```

Here the squared norm $\tfrac{1}{2} \| \mathbf{A} \mathbf{x} - \mathbf{y} \|^2$ is a *smooth* function $f$ whereas the $l_1$-norm is a *nonsmooth* function $g$. This problem can be solved with only few lines of code:

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

!!! note

    The function `ls` is a short hand notation for `0.5*norm(...)^2`, namely a least squares term.


It is possible to access to the solution by typing `~x`.
By default variables are initialized by `Array`s of zeros.
Different initializations can be set during construction `x = Variable( [1.; 0.; ...] )` or by assignment `~x .= [1.; 0.; ...]`.

## Constrained optimization

Constrained optimization is also encompassed by the [Standard problem formulation](@ref): for a nonempty set $\mathcal{S}$ the constraint of

```math
\begin{align*}
\underset{ \mathbf{x} }{\text{minimize}} \ &  f(\mathbf{x}) \\
\text{subject to} \ & \mathbf{x} \in \mathcal{S}
\end{align*}
```

can be converted into an *indicator function*

```math
g(\mathbf{x}) = \delta_{\mathcal{S}} (\mathbf{x}) =  \begin{cases}
    0       & \text{if} \ \mathbf{x} \in \mathcal{S},\\
    +\infty & \text{otherwise}.
    \end{cases}
```

Constraints are treated as *nonsmooth functions*.
This conversion is automatically performed by StructuredOptimization.jl.
For example, the non-negative deconvolution problem:

```math
\begin{align*}
\underset{ \mathbf{x} }{\text{minimize}} \ &  \tfrac{1}{2} \| \mathbf{x} * \mathbf{h} - \mathbf{y} \|^2 \\
\text{subject to} \ & \mathbf{x} \geq 0
\end{align*}
```

where $*$ stands for convolution and $\mathbf{h}$ contains the taps of a finite impulse response filter,
can be solved using the following lines of code:

```julia
julia> n = 10;                        # define problem size

julia> x = Variable(n);               # define variable

julia> h, y = randn(n), randn(2*n-1); # random filter taps and output

julia> @minimize ls(conv(x,h)-y) st x >= 0.

```

!!! note

    The convolution mapping was applied to the variable `x` using `conv`.
    StructuredOptimization.jl provides a set of functions that can be
    used to apply specific operators to variables and create mathematical
    expression. The available functions can be found in [Mappings](@ref).
    In general it is more convenient to use these functions instead of matrices,
    as these functions apply efficient algorithms for the forward and adjoint
    mappings leading to *matrix free optimization*.

## Using multiple variables

It is possible to use multiple variables which are allowed to be matrices or even tensors. For example a non-negative matrix factorization problem:

```math
\begin{align*}
\underset{ \mathbf{X}_1, \mathbf{X}_2  }{\text{minimize}} \ &  \tfrac{1}{2} \| \mathbf{X}_1 \mathbf{X}_2 - \mathbf{Y} \| \\
\text{subject to} \ & \mathbf{X}_1 \geq 0,  \ \mathbf{X}_2 \geq 0,
\end{align*}
```

can be solved using the following code:

```julia
# matrix variables initialized with random coefficients
julia> X1, X2 = Variable(rand(n,l)), Variable(rand(l,m));

julia> Y = rand(n,m);

julia> @minimize ls(X1*X2-Y) st X1 >= 0., X2 >= 0.

```

## Limitations

Currently StructuredOptimization.jl supports only *proximal gradient algorithms* (i.e., *forward-backward splitting* base), which require specific properties of the nonsmooth functions and constraint to be applicable. In particular, the nonsmooth functions must have an *efficiently computable proximal mapping*.

If we express the nonsmooth function $g$ as the composition of
a function $\tilde{g}$ with a linear operator $A$:

```math
g(\mathbf{x}) =
\tilde{g}(A \mathbf{x})
```

then the proximal mapping of $g$ is efficiently computable if either of the following hold:

1. Operator $A$ is a *tight frame*, namely it satisfies $A A^* = \mu Id$, where $\mu \geq 0$, $A^*$ is the adjoint of $A$, and $Id$ is the identity operator.

2. Function $g$ is a *separable sum* $g(\mathbf{x}) = \sum_j h_j (B_j \mathbf{x}_j)$, where $\mathbf{x}_j$ are non-overlapping slices of $\mathbf{x}$, and $B_j$ are tight frames.

Let us analyze these rules with a series of examples.
The LASSO example above satisfy the first rule:

```julia
julia> @minimize ls( A*x - y ) + λ*norm(x, 1)
```

since the nonsmooth function $\lambda \| \cdot \|_1$ is not composed with any operator (or equivalently is composed with $Id$ which is a tight frame).
Also the following problem would be accepted by StructuredOptimization.jl:

```julia
julia> @minimize ls( A*x - y ) + λ*norm(dct(x), 1)
```

since the discrete cosine transform (DCT) is orthogonal and is therefore a tight frame. On the other hand, the following problem

```julia
julia> @minimize ls( A*x - y ) + λ*norm(x, 1) st x >= 1.0
```

cannot be solved through proximal gradient algorithms, since the second rule would be violated.
Here the constraint would be converted into an indicator function and the nonsmooth function $g$ can be written as the sum:

```math
g(\mathbf{x}) =\lambda \| \mathbf{x} \|_1 + \delta_{\mathcal{S}} (\mathbf{x})
```

which is not separable. On the other hand this problem would be accepted:

```julia
julia> @minimize ls( A*x - y ) + λ*norm(x[1:div(n,2)], 1) st x[div(n,2)+1:n] >= 1.0
```

as not the optimization variables $\mathbf{x}$ are partitioned into non-overlapping groups.

!!! note

    When the problem is not accepted it might be still possible to solve it: see [Smoothing](@ref) and [Duality](@ref).
