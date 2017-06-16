export solve

abstract type Solver end
abstract type ForwardBackwardSolver <: Solver end

include("solvers/pg.jl")
include("solvers/zerofpr.jl")
include("solvers/utils.jl")

default_slv = FPG #TODO change this to ZeroFPR once ready

# To print out solver objects

function Base.show(io::IO, slv::ForwardBackwardSolver)
	println(io, fun_name(slv) )
	println(io, "iterations : $(slv.it) / $(slv.maxit)")
	println(io, "fpr        : $(deepmaxabs(slv.normfpr)/slv.gamma)")
	println(io, "cost       : $(slv.cost)")
	println(io, "gamma      : $(slv.gamma)")
	println(io, "time       : $(slv.time)")
	println(io, "prox   eval: $(slv.cnt_prox)")
	print(  io, "matvec eval: $(slv.cnt_matvec)")
end

