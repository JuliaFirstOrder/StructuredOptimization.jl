# AbstractOperators binding
# special cases
import Base: *, reshape

"""
`reshape(ex::AbstractExpression, dims...)`

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
	if typeof(displacement(A)) <: Number
		d = displacement(A)
	else
		d = reshape(displacement(A), dims...)
	end
	return Expression{length(A.x)}(A.x,op,d)
end
#Reshape

imported = [:getindex :GetIndex;
            :fft      :DFT;
            :rfft     :RDFT;
            :irfft    :IRDFT;
            :ifft     :IDFT;
            :dct      :DCT;
            :idct     :IDCT;
            :conv     :Conv;
            :xcorr    :Xcorr;
            :filt     :Filt;
           ]

exported = [:finitediff :FiniteDiff;
            :variation  :Variation;
            :mimofilt   :MIMOFilt;
            :zeropad    :ZeroPad;
            :sigmoid    :Sigmoid;
            :σ          :Sigmoid; #alias
           ]

#importing functions from Base
for f in  imported[:,1]
	@eval begin
		import Base: $f
	end
end
#exporting functions
for f in  exported[:,1]
	@eval begin
		export $f
	end
end

fun = [imported; exported]
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
`getindex(ex::AbstractExpression, dims...)`

Slices the codomain of `ex`.

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> x[1:4]

julia> getindex(randn(15,10)*x,1:10)

```
"""
getindex

"""
`fft(ex::AbstractExpression)`

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
`rfft(ex::AbstractExpression, [, dims] )`

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
`ifft(ex::AbstractExpression)`

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
`irfft(ex::AbstractExpression, d, [, dims] )`

Applies the inverse Fourier transform exploiting the conjugate symmetry. 
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
`dct(ex::AbstractExpression)`

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
`idct(ex::AbstractExpression)`

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
`conv(ex::AbstractExpression, h::AbstractVector)`

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
`xcorr(ex::AbstractExpression, h::AbstractVector)`

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
`filt(ex::AbstractExpression, b::AbstractVector, [a::AbstractVector])`

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
`finitediff(ex::AbstractExpression, dir = 1)`

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
`variation(ex::AbstractExpression)`

Performs the variation using forward finite differences. 
Specifically if the codomain of `ex` is  `ℝ^(n, m)` the codomain of `variation(ex)` will 
consist of `ℝ^(n*m,2)`, having in the `i`th column the vectorized gradient over the `i`th direction. 

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
`mimofilt(ex::AbstractExpression, B::Vector{AbstractVector}, [A::Vector{AbstractVector}])`

Multiple-input multiple-output filter: like `filt` but allows for multipule inputs. 

See documentation of `AbstractOperator.MIMOFilt`.  
"""
mimofilt   

"""
`zeropad(ex::AbstractExpression, zp::Tuple)`

Performs zeropadding:
```math
zp(\\mathbf{x}) = [\\mathbf{x}, 0, \\dots, 0  ]^T
```
where `zp` is the number of nonzero elements that are added.

See documentation of `AbstractOperator.ZeroPad`.  
"""
zeropad    

"""
`sigmoid(ex::AbstractExpression, γ = 1.0)`
`σ(ex::AbstractExpression, γ = 1.0)`

Sigmoid function:
```math
\\sigma(\\mathbf{x}) = \\frac{1}{1+e^{-\\gamma \\mathbf{x} } }
```

See documentation of `AbstractOperator.Sigmoid`.  
"""
sigmoid    
σ          
