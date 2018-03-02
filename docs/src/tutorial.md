# Quick Tutorial

## Standard problem formulation

Currently with `StructuredOptimization.jl` you can solve problems of the form

```math
\underset{ \mathbf{x} }{\text{minimize}} \ f(\mathbf{x}) + g(\mathbf{x}),
```

where $f$ is a smooth function while $g$ is possibly nonsmooth.

## Unconstraint optimization

The LASSO problem is popular example of this class of problems: 

```math
\underset{ \mathbf{x} }{\text{minimize}} \ \tfrac{1}{2} \| \mathbf{A} \mathbf{x} - \mathbf{y} \|^2+  \| \mathbf{x} \|_1.
```

Here the squared norm $\tfrac{1}{2} \| \mathbf{A} \mathbf{x} - \mathbf{y} \|^2$ is a _smooth_ function while the $l_1$-norm is a _nonsmooth_ function.

This can be solved using `StructuredOptimization.jl` using only few lines of code:

```julia
julia> using StructuredOptimization

julia> n, m = 100, 10;                # define problem size 

julia> A, y = randn(m,n), randn(m);   # random problem data

julia> x = Variable(n);               # initialize optimization variable

julia> λ = 1e-2*norm(A'*y,Inf);       # define λ    

julia> @minimize ls( A*x - y ) + λ*norm(x, 1); # minimize problem

```

!!! note 

    The function `ls` is a short hand notation for `0.5*norm(...)^2`, namely a least squares term.


It is possible to access to the solution by typing `~x`. 

By default variables are initialized by `Array`s of zeros. 

It is possible to set different initializations during construction `x = Variable( [1.; 0.; ...] )` or by assignement `~x .= [1.; 0.; ...]`.

## Constraint optimization

Constraint optimization is also ecompassed by [Standard problem formulation](@ref): 

for a nonempty set $\mathcal{S}$ the constraint of 

```math
\begin{align*}
\underset{ \mathbf{x} }{\text{minimize}} \ &  f(\mathbf{x}) \\
\text{subject to} \ & \mathbf{x} \in \mathcal{S}
\end{align*}
```

can be converted into an indicator function

```math
g(\mathbf{x}) =  \begin{cases}
    0       & \text{if} \ \mathbf{x} \in \mathcal{S},\\
    +\infty & \text{otherwise},
    \end{cases}
```

to obtain the standard form. Constraints are treated as _nonsmooth functions_.

This conversion is automatically performed by `StructuredOptimization.jl`.

For example, the non-negative deconvolution problem:

```math
\begin{align*}
\underset{ \mathbf{x} }{\text{minimize}} \ &  \tfrac{1}{2} \| \mathbf{x} * \mathbf{h} - \mathbf{y} \| \\
\text{subject to} \ & \mathbf{x} \geq 0
\end{align*}
```

where $*$ stands fof convoluton and $\mathbf{h}$ contains the taps of a finite impluse response.

This problem be solved using the following line of code:

```julia
julia> n = 10;

julia> x = Variable(n);               # define variable

julia> h, y = randn(n), randn(2*n-1); # random filter taps and output

julia> @minimize ls(conv(x,h)-y) st x >= 0.

```

!!! note 

    The convolution mapping was applied to the variable `x` using `conv`. 

    `StructuredOptimization.jl` provides a set of functions that can be used to apply 
    specific operators to variables and create mathematical expression. 
    
    The available functions can be found in [Operators](@ref).

## Using multiple variables

It is possible to use multiple variables which are allowed to be matrices or even tensors. 

For example a non-negative matrix factorization problem:

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

**TODO simplify this**

Currently `StructuredOptimization.jl` supports only Proximal Gradient (aka Forward Backward) algorithms, which require certain properties of the nonsmooth functions and costraint.

In the general case a nonsmooth function of $M$ variables composed by $G$ terms can be written as: 
```math
g(\mathbf{x}_1,\dots,\mathbf{x}_M) =
\sum_{i = 0}^G g_i \left(\sum_{j = 1}^{M}
A_{i,j} \mathbf{x}_j \right).
```
where the functions $g_i$ are nonsmooth functions (or indicator functions resulting from constraints) and $A_{i,j}$ linear operators.

The problem can be solved when $g$ satisfies the following conditions:

1. for all $i\in \{1,\ldots,G \}$ and $j\in\{1,\ldots,M \}$, mapping $A_{i,j}$ satisfies $A_{i,j}^* A_{i,j} = \mu_{i,j} I$, where $\mu_{i,j} \geq 0$, $A^*$ is the adjoint of $A$ and $\mathcal{I}$ is the identity operator.

2. for all $j \in \{1,\dots,M \}$, the cardinality of $\{i | A_{i,j} \neq 0 \} = 1$. 

Let us analyze these rules with a series of examples. 

The previous example was satisfing the rules:
```julia
@minimize ls(X1*X2-Y) st X1 >= 0., X2 >= 0.
```
Here there are two constraints each one containing only one variable and 







