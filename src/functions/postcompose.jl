is_smooth(f::Postcompose) = is_smooth(f.f)
is_quadratic(f::Postcompose) = is_quadratic(f.f)
is_generalized_quadratic(f::Postcompose) = is_generalized_quadratic(f.f)
is_strongly_convex(f::Postcompose) = is_strongly_convex(f.f)
