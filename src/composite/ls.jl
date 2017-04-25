export ls

ls(x::Variable) = ls(eye(x))
ls(A::AbstractAffineTerm) = CompositeFunction(variable(A), [SqrNormL2(1.)], [A])
