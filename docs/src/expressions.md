# Expressions

With `StructuredOptimization.jl` you can easily create mathematical expressions.
Firstly, [Variables](@ref) must be defined: various [Mappings](@ref) can then
be applied following the application of [Functions and constraints](@ref) to
create the `Term`s  that define the optimization problem.

## Variables

### Creating Variables

```@docs
Variable
```
!!! note

    `StructuredOptimization.jl` supports complex variables. It is possible to create them by specifying the type
    `Variable(Complex{Float64}, 10)` or by initializing them with a complex array `Variable(randn(10)+im*randn(10))`.

### Utilities

```@docs
~
size
eltype
```

## Summing expressions

```@docs
+
```

## Multiplying expressions

```@docs
*
```

## Mappings

As shown in the [Quick tutorial guide](@ref) it is possible to apply different mappings to the variables
using a simple syntax.

Alternatively, as shown in [Multiplying expressions](@ref), it is possible to define the mappings using
[`AbstractOperators.jl`](https://github.com/kul-forbes/ProximalAlgorithms.jl) and to apply them
to the variable (or expression) through multiplication.

### Basic mappings
```@docs
getindex
reshape
```

### DSP mappings
```@docs
fft
ifft
rfft
irfft
dct
idct
conv
xcorr
filt
mimofilt
zeropad
```

### Finite differences mappings
```@docs
finitediff
variation
```

### Nonlinear mappings
```@docs
sin
cos
atan
tanh
exp
pow
sigmoid
```

## Utilities

It is possible to access the variables, mappings and displacement of an expression.
Notice that these commands work also for the `Term`s described in [Functions and constraints](@ref).

```@docs
variables
operator
affine
AbstractOperators.displacement
```
