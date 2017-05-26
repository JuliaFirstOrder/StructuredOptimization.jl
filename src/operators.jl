export LinearOperator

abstract type LinearOperator end

import Base: A_mul_B!, Ac_mul_B!, size, ndims, transpose, *, +, -

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
# include("operators/LBFGS.jl")
# include("operators/BlkDiagLBFGS.jl")
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

is_null(L::LinearOperator) = false
is_eye(L::LinearOperator) = false
is_diagonal(L::LinearOperator) = false
is_gram_diagonal(L::LinearOperator) = is_diagonal(L)
is_invertible(L::LinearOperator) = false
is_full_row_rank(L::LinearOperator) = false
is_full_column_rank(L::LinearOperator) = false

function Base.show(io::IO, L::LinearOperator)
  print(io, typeof(L))
end

# Shorthands

function (*){T <: Union{AbstractArray, Tuple}}(L::LinearOperator, b::T)
  y = deepzeros(codomainType(L), size(L, 1))
	A_mul_B!(y, L, b)
  return y
end

(-){T <: LinearOperator}(L::T) = Scale(-1.0, L)

(+){T <: LinearOperator}(L::T) = L

transpose{T <: LinearOperator}(L::T) = Transpose(L)
