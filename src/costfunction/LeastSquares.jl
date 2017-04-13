export ls

ls(x::Variable)       = ls(eye(x))
ls(A::AffineOperator) = CostFunction(variable(A), [LeastSquares(1.)], [A])

