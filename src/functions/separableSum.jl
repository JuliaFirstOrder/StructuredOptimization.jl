is_smooth(f::SeparableSum) = all(is_smooth.(f.fs))
is_quadratic(f::SeparableSum) = all(is_quadratic.(f.fs))
is_strongly_convex(f::SeparableSum) = all(is_strongly_convex.(f.fs))

function gradient!(grad::AbstractArray, fs::AbstractArray, x::AbstractArray)
  val = 0.0
  for k in eachindex fs
    val += gradient!(grad[k], fs[k], x[k])
  end
  return val
end

gradient!(grad::AbstractArray, f::SeparableSum, x::AbstractArray) =
  gradient!(grad, f.fs, x)
