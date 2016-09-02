function print_status(it, gamma, fpr, cost, verbose)
  if verbose > 0 && it == 1
    @printf("%6s | %10s | %10s | %14s\n", "it", "gamma", "fpr", "cost")
    @printf("-------|------------|------------|---------------\n")
  end
  if verbose == 2 || (verbose == 1 && (it == 1 || it%100 == 0))
    @printf("%6d | %7.4e | %7.4e | %11.8e\n", it, gamma, fpr, cost)
  end
end

function halt(tol::Float64, gamma::Float64, fpr0::Float64, fpr::Float64, fun_prev::Float64, fun::Float64)
	conv_fpr = fpr <= (1+fpr0)*tol
	conv_fun = abs(fun-fun_prev) <= (1+abs(fun))*tol
	return conv_fpr && conv_fun
end

abstract SolverInfo

immutable BasicInfo <: SolverInfo
  iterations::Int
  gamma::Float64
  fpr::Float64
  cost::Float64
  time::Float64
end

function Base.show(io::IO, info::BasicInfo)
  println(io, "iterations : $(info.iterations)")
  println(io, "gamma      : $(info.gamma)")
  println(io, "fpr        : $(info.fpr)")
  println(io, "cost       : $(info.cost)")
  print(  io, "time       : $(info.time)")
end
