import Base: rank

rank(x::Variable, args...) = rank(eye(x), args...)

rank(A::AffineOperator) = CostFunction(variable(A), [IndBallRank(0)], [A])

<=(T::IndBallRank,   r::Integer) = IndBallRank(r) 

