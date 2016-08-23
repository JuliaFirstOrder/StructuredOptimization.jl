# Definitions useful to work with tuples of variable blocks.
# First we extend sum, subtraction, to (omogeneous) tuples.

import Base.+, Base.-

function +(a::Tuple, b::Tuple)
   N = length(a)
   if length(b) != N error("operands must have the same length") end
   return map(i -> a[i]+b[i], (1:N...))
end

function -(a::Tuple, b::Tuple)
   N = length(a)
   if length(b) != N error("operands must have the same length") end
   return map(i -> a[i]-b[i], (1:N...))
end

# Then we would need a way to stack linear operators horizontally/vertically
# pretty much like they do in LinearOperators, but they only handle operators
# from/to spaces of Array{T,1}.

# Next we define the separable sum of proximable functions.
# This is basically the sum of a tuple of functions, each one applied to the
# correspondent element of a tuple of variables (of appropriate dimension).

immutable funSum <: ProximableFunction
  N::Int
  fs::Tuple
  funSum(fs...) =
    new(length(fs), fs)
end

function call(f::funSum, xs...)
  vs = map(i -> f.fs[i](xs[i]), 1:f.N)
  return sum(vs)
end

function prox(f::funSum, gamma::Float64, xs...)
  res = map(i -> prox(f.fs[i], gamma, xs[i]), (1:f.N...))
  return map(p -> p[1], res), sum(map(p -> p[2], res))
end

function Base.show(io::IO, f::funSum)
  for i = 1:f.N-1
    println(io, f.fs[i])
    println(io, "---")
  end
  print(io, f.fs[f.N])
end
