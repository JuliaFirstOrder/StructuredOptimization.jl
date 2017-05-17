export LinearOperator

abstract type LinearOperator end

import Base: A_mul_B!, Ac_mul_B!, size, ndims
  # *,
  # .*,
  # transpose,
  # inv,
  # size,
  # ndims,
  # +,
  # -,
  # ==

# +(L::LinearOperator) = L

# function *(L::LinearOperator,b::AbstractArray)
# 	y = zeros(codomainType(L),size(L,1))
# 	A_mul_B!(y,L,b)
# 	return y
# end

# Basic operators

include("operators/MyOperator.jl")
include("operators/Zeros.jl")
include("operators/Eye.jl")
include("operators/DiagOp.jl")
include("operators/GetIndex.jl")
include("operators/MatrixOp.jl")
include("operators/DFT.jl")
include("operators/DCT.jl")
include("operators/FiniteDiff.jl")
include("operators/Variation.jl")
include("operators/Conv.jl")
include("operators/Filt.jl")
include("operators/MIMOFilt.jl")
include("operators/ZeroPad.jl")
include("operators/Xcorr.jl")
include("operators/LBFGS.jl")
include("operators/BlkDiagLBFGS.jl")
include("operators/utils.jl")

# Calculus rules

include("operators/calculus/DCAT.jl")
include("operators/calculus/HCAT.jl")
include("operators/calculus/VCAT.jl")
include("operators/calculus/Compose.jl")
include("operators/calculus/Reshape.jl")
include("operators/calculus/Scale.jl")
include("operators/calculus/Sum.jl")
include("operators/calculus/Transpose.jl")

size(L::LinearOperator, i::Int) = size(L)[i]
ndims(L::LinearOperator) = length(size(L,1)), length(size(L,2))
ndims(L::LinearOperator, i::Int) = ndims(L)[i]

#usually domain is preserved, if not one has to redefine these functions
  domainType(L::LinearOperator) = L.domainType
codomainType(L::LinearOperator) = L.domainType

is_diagonal(L::LinearOperator) = false
is_gram_diagonal(L::LinearOperator) = is_diagonal(L::LinearOperator)
is_invertible(L::LinearOperator) = false
is_full_row_rank(L::LinearOperator) = false
is_full_column_rank(L::LinearOperator) = false

function Base.show(io::IO, L::LinearOperator)
  println(io, "description : ", fun_name(L))
  print(  io, "type        : ", fun_type(L))
end

fun_type(  L) = fun_domain(L)*" → "*fun_codomain(L)

fun_domain(L::LinearOperator)   =   domainType(L) <: Complex ? "ℂ^$(size(L,2))" : "ℝ^$(size(L,2))"
fun_codomain(L::LinearOperator) = codomainType(L) <: Complex ? "ℂ^$(size(L,1))" : "ℝ^$(size(L,1))"
