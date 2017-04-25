export ls

ls(x::Variable)       = ls(eye(x))
ls(A::AbstractAffineTerm) = CompositeFunction(variable(A), [LeastSquares(1.)], [A])

