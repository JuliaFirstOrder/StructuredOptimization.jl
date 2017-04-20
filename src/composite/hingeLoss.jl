export hingeloss

hingeloss(x::Variable, args...) = hingeloss(eye(x), args...)

hingeloss{R <: Real}(A::AbstractAffineTerm, b::Array{R,1}) =
CompositeFunction(variable(A), [HingeLoss(b, 1.0)], [A])
