is_smooth(f::SeparableSum) = all(is_smooth.(f.fs))
is_quadratic(f::SeparableSum) = all(is_quadratic.(f.fs))
is_strongly_convex(f::SeparableSum) = all(is_strongly_convex.(f.fs))

function gradient!{T <: Tuple}(grad::T, fs::Tuple, x::T)
  val = 0.0
  for k in eachindex(fs)
    val += gradient!(grad[k], fs[k], x[k])
  end
  return val
end

gradient!{T <: Tuple}(grad::T, f::SeparableSum, x::T) =
  gradient!(grad, f.fs, x)
