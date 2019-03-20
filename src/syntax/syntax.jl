abstract type AbstractExpression end

include("variable.jl")
include("expressions/expression.jl")
include("terms/term.jl")
include("problem.jl")

const TermOrExpr =  Union{Term,AbstractExpression}
