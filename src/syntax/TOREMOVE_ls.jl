export ls

# ls(x::Variable) = ls(LinearExpression(Eye()))
# ls(A::AbstractAffineExpression) = Term(variable(A), [SqrNormL2(1.)], [A])

ls(ex) = Term(SqrNormL2(), ex)
