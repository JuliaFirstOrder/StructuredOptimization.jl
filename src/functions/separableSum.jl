immutable SeparableSum <: ExtendedRealValuedFunction
	fs::ExtendedRealValuedFunction
end

function (fs::AbstractArray{ExtendedRealValuedFunction})(x::AbstractArray)
  val = 0.0
  for k in eachindex fs
    val += fs[k](x[k])
  end
  return val
end

(f::SeparableSum)(x::AbstractArray) = f.fs(x)

function gradient!(grad::AbstractArray, fs::AbstractArray{ExtendedRealValuedFunction}, x::AbstractArray)
  val = 0.0
  for k in eachindex fs
    val += gradient!(grad[k], fs[k], x[k])
  end
  return val
end
