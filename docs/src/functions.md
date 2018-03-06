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
logbarrier
```

## Nonsmooth functions

```@docs
norm
maximum
sumpositive
hingeloss
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

Sometimes the optimization problem might involve only non-smooth terms which do not lead to efficient proximal mappings. It is possible to *smooth* this terms by means of the *Moreau envelope*.

```@docs
smooth
```

## Duality

In some cases it is more convenient to solve the *dual problem* instead of the primal problem. 

It is possible to convert the primal problem into its dual form by means of the *convex conjugate*. 

See the Total Variation demo for an example of such procedure.

```@docs
conj
```

