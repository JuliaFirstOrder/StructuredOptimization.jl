export hingeloss

hingeloss(x::Variable, args...) = hingeloss(eye(x), args...)

hingeloss{R <: Real}(A::AffineOperator, b::Array{R,1}) = 
CostFunction(variable(A), [HingeLoss(b, 1.0)], [A])

