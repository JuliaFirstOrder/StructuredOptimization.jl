export hingeloss

hingeloss(x::Variable, args...) = hingeloss(eye(x), args...)

hingeloss{R <: Real}(A::AbstractAffineExpression, b::Array{R,1}) =
Term(variable(A), [HingeLoss(b, 1.0)], [A])
