function gradient!{T <: RealOrComplex}(y::AbstractArray{T}, g::PrecomposeDiagonal, x::AbstractArray{T})
  z = g.a .* x .+ g.b
  v = gradient!(y, g.f, z)
  y .*= g.a
  return v
end
