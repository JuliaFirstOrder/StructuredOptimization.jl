rank(x::Variable, args...) = rank(eye(x), args...)

# Dirty trick: the "rank" function only makes sense in constraints such as
#   rank(X) <= r,
# therefore here the parameter (1) doesn't really have a role.
# We should probably fix this: it allows weird things in expressing problems.
# Maybe we should have Rank <: ProximableFunction (with no prox! nor gradient!
# defined), that gives IndBallRank when combined with <=.
rank(A::AbstractAffineTerm) = CompositeFunction(variable(A), [IndBallRank(1)], [A])
