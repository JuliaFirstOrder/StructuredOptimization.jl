is_smooth(f::Conjugate) = is_strongly_convex(f.f)
is_quadratic(f::Conjugate) = is_strongly_convex(f.f) && is_quadratic(f.f)
is_generalized_quadratic(f::Conjugate) = is_quadratic(f.f)
is_strongly_convex(f::Conjugate) = is_smooth(f.f)
