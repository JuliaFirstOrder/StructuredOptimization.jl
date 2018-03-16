# Functions and constraints

Once an expression is created it is possible to create the `Term`s defining the optimization problem.
These can consists of either [Smooth functions](@ref),  [Nonsmooth functions](@ref), [Inequality constraints](@ref)
or [Equality constraints](@ref).

## Smooth functions

```@docs
ls
huberloss
sqrhingeloss
crossentropy
logisticloss
dot
```

## Nonsmooth functions

```@docs
norm
maximum
sumpositive
hingeloss
logbarrier
```

## Inequality constraints

```@docs
<=
```

## Equality constraints

```@docs
==
```

## Smoothing

Sometimes the optimization problem might involve non-smooth terms which
do not have efficiently computable proximal mappings.
It is possible to *smoothen* these terms by means of the *Moreau envelope*.

```@docs
smooth
```

## Duality

In some cases it is more convenient to solve the *dual problem* instead
of the primal problem. It is possible to convert a problem into its dual
by means of the *convex conjugate*.

See the [Total Variation demo](https://github.com/kul-forbes/StructuredOptimization.jl/blob/master/demos/TotalVariationDenoising.ipynb) for an example of such procedure.

```@docs
conj
```
