# AbstractOperators binding
# special cases
import Base: *, reshape

"""
    reshape(x::AbstractExpression, dims...)

Reshapes the codomain of the expression.

# Example

```julia
julia> A,b = randn(10,3), randn(10);

julia> reshape(A*x-b,2,5)

```
"""
function reshape(a::AbstractExpression, dims...)
  A = convert(Expression,a)
  op = Reshape(A.L, dims...)
  return Expression{length(A.x)}(A.x,op)
end
#Reshape

imported = [
            :getindex :GetIndex;
            :exp      :Exp;
            :cos      :Cos;
            :sin      :Sin;
            :atan     :Atan;
            :tanh     :Tanh;
           ]

importedFFTW = [
                :fft      :(AbstractOperators.DFT);
                :rfft     :RDFT;
                :irfft    :IRDFT;
                :ifft     :IDFT;
                :dct      :DCT;
                :idct     :IDCT;
               ]

importedDSP = [
               :conv     :Conv;
               :xcorr    :Xcorr;
               :filt     :Filt;
              ]

exported = [
            :finitediff :FiniteDiff;
            :variation  :Variation;
            :mimofilt   :MIMOFilt;
            :zeropad    :ZeroPad;
            :sigmoid    :Sigmoid;
            :σ          :Sigmoid; #alias
            :pow        :Pow; #alias
           ]

#importing functions from Base
for f in  imported[:,1]
  @eval begin
    import Base: $f
  end
end
#importing functions from FFTW
for f in  importedFFTW[:,1]
  @eval begin
    import FFTW: $f
    export $f
  end
end
#importing functions from DSP
for f in  importedDSP[:,1]
  @eval begin
    import DSP: $f
    export $f
  end
end
#exporting functions
for f in  exported[:,1]
  @eval begin
    export $f
  end
end

fun = [imported; importedFFTW; importedDSP; exported]
for i = 1:size(fun,1)
  f,fAbsOp = fun[i,1],fun[i,2]
  @eval begin
    function $f(a::AbstractExpression, args...)
      A = convert(Expression,a)
      op = $fAbsOp(codomainType(operator(A)),size(operator(A),1), args...)
      return op*A
    end
  end
end

## docs

"""
    getindex(x::AbstractExpression, dims...)

Slices the codomain of `ex`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> getindex(x,1:4)

julia> X = Variable(10,4)
Variable(Float64, (10, 4))

julia> X[:,1:3]

```
"""
getindex

"""
    fft(x::AbstractExpression)

Applies the Fourier transform.
See documentation for `Base.fft` and `AbstractOperator.DFT`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = fft(x)

julia> operator(ex)
ℱ  ℝ^10 -> ℂ^10

```
"""
fft

"""
    rfft(x::AbstractExpression, [, dims] )

Applies the Fourier transform exploiting the conjugate symmetry.
See documentation for `Base.rfft` and `AbstractOperator.RDFT`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = rfft(x)

julia> operator(ex)
ℱ  ℝ^10 -> ℂ^6

```
"""
rfft

"""
    ifft(x::AbstractExpression)

Applies the inverse Fourier transform.
See documentation for `Base.ifft` and `AbstractOperator.IDFT`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = ifft(x)

julia> operator(ex)
ℱ⁻¹  ℝ^10 -> ℂ^10

```
"""
ifft



"""
    irfft(x::AbstractExpression, d, [, dims] )

Applies the inverse Fourier transform exploiting the conjugate symmetry.

`d` must indicate the length of the real codomain.

See documentation for `Base.irfft` and `AbstractOperator.IRDFT`.

# Example

```julia
julia> x = Variable(Complex{Float64}, 10)
Variable(Float64, (10,))

julia> ex = irfft(x,19)

julia> operator(ex)
ℱ⁻¹  ℂ^10 -> ℝ^19

```
"""
irfft

"""
    dct(x::AbstractExpression)

Applies the discrete cosine transform.
See documentation for `Base.dct` and `AbstractOperator.DCT`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = dct(x)

julia> operator(ex)
ℱc  ℝ^10 -> ℝ^10

```
"""
dct

"""
    idct(x::AbstractExpression)

Applies the inverse discrete cosine transform.
See documentation for `Base.idct` and `AbstractOperator.IDCT`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = idct(x)

julia> operator(ex)
ℱc⁻¹  ℝ^10 -> ℝ^10

```
"""
idct

"""
    conv(x::AbstractExpression, h::AbstractVector)

Perform discrete convolution with `h`.
See documentation for `Base.conv` and `AbstractOperator.Conv`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = conv(x, randn(5))

julia> operator(ex)
★  ℝ^10 -> ℝ^14

```
"""
conv

"""
    xcorr(x::AbstractExpression, h::AbstractVector)

Performs the cross correlation of the codomain of `ex` with `h`.
See documentation for `Base.xcorr` and `AbstractOperator.Xcorr`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = xcorr(x, randn(5))

julia> operator(ex)
◎  ℝ^10 -> ℝ^19

```
"""
xcorr

"""
    filt(x::AbstractExpression, b::AbstractVector, [a::AbstractVector])

Filter with the finite impulse response (FIR) `b`.
Alternatively infinite impulse responses filters (IIR) can be used as well by specifying the coefficients `a`.

See documentation for `Base.filt` and `AbstractOperator.Filt`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = filt(x, [1.;0.;0.],[1.;0.;1.])

julia> operator(ex)
IIR  ℝ^10 -> ℝ^10

julia> ex = filt(x, [1.;0.;0.])

julia> operator(ex)
FIR  ℝ^10 -> ℝ^10

```
"""
filt

"""
    finitediff(X::AbstractExpression, dir = 1)

Performs the discretized gradient over the specified direction `dir` obtained using forward finite differences.

# Example

```julia
julia> X = Variable(10,10)
Variable(Float64, (10,10))

julia> ex = finitediff(X)

julia> operator(ex)
δx  ℝ^(10, 10) -> ℝ^(9, 10)

julia> ex = finitediff(X,2)

julia> operator(ex)
δy  ℝ^(10, 10) -> ℝ^(10, 9)

```
"""
finitediff

"""
    variation(X::AbstractExpression)

Returns the variation of `X` using forward finite differences.
Specifically if `X` is in `ℝ^(n, m)` the codomain of `variation(X)` will
consist of `ℝ^(n*m,2)`, having in the ``i``-th column the vectorized gradient over the ``i``-th direction.

See documentation of `AbstractOperator.Variation`.

# Example

```julia
julia> X = Variable(4,5,6)
Variable(Float64, (4,5,6))

julia> ex = variation(X)

julia> operator(ex)
Ʋ  ℝ^(4, 5, 6) -> ℝ^(120, 3)

```
"""
variation

"""
    mimofilt(X::AbstractExpression, B::Vector{AbstractVector}, [A::Vector{AbstractVector}])

Multiple-input multiple-output filter: like `filt` but allows for multipule inputs.

```math
\\mathbf{y}_i = \\sum_{j = 1}^{M} \\mathbf{h}_{i,j} * \\mathbf{x}_j
```
where ``\\mathbf{y}_i`` and ``\\mathbf{x}_j`` are the ``i``-th and ``j``-th columns of the output `Y` and input `X` matrices respectively.

The filters ``\\mathbf{h}_{i,j}`` can be represented either by providing coefficients `B` and `A` (IIR) or `B` alone (FIR). These coefficients must be given in a `Vector` of `Vector`s.

For example for a `3` by `2` MIMO system (i.e. `size(X,2) == 3` inputs and `size(Y,2) == 2` outputs) `B` must be:

`B = [b11, b12, b13, b21, b22, b23]`

where `bij` are vector containing the filter coeffients of `h_{i,j}`.

See documentation of `AbstractOperator.MIMOFilt`.
"""
mimofilt

"""
    zeropad(x::AbstractExpression, zp::Tuple)

Performs zeropadding:
```math
zp(\\mathbf{x}) = [\\mathbf{x}, 0, \\dots, 0  ]^T
```
where `zp` is a tuple containing the number of nonzero elements that are added in each dimension.

See documentation of `AbstractOperator.ZeroPad`.
"""
zeropad

"""
    sigmoid(x::AbstractExpression, γ = 1.0)
`σ(x::AbstractExpression, γ = 1.0)`

Sigmoid function:
```math
\\sigma(\\mathbf{x}) = \\frac{1}{1+e^{-\\gamma \\mathbf{x} } }
```

See documentation of `AbstractOperator.Sigmoid`.
"""
sigmoid
σ

"""
    exp(x::AbstractExpression)

Exponential function:
```math
e^{ \\mathbf{x}  }
```

See documentation of `AbstractOperator.Exp`.
"""
exp

"""
    sin(x::AbstractExpression)

Sine function:
```math
\\sin( \\mathbf{x}  )
```

See documentation of `AbstractOperator.Sin`.
"""
sin

"""
    cos(x::AbstractExpression)

Cosine function:
```math
\\cos( \\mathbf{x}  )
```

See documentation of `AbstractOperator.Cos`.
"""
cos

"""
    atan(x::AbstractExpression)

Inverse tangent function:
```math
\\tan^{-1}( \\mathbf{x}  )
```

See documentation of `AbstractOperator.Atan`.
"""
atan

"""
    tanh(x::AbstractExpression)

Hyperbolic tangent function:
```math
\\tanh ( \\mathbf{x}  )
```

See documentation of `AbstractOperator.Tanh`.
"""
tanh

"""
    pow(x::AbstractExpression, n)

Elementwise `n`-th power of `x`:
```math
x_i^{n} \\ \\forall \\  i = 0,1, \\dots
```

See documentation of `AbstractOperator.Pow`.
"""
pow
