function print_status(slv::ForwardBackwardSolver)
  if slv.verbose > 0 && slv.it == 1
    @printf("%6s | %10s | %10s | %14s\n", "it", "gamma", "fpr", "cost")
    @printf("-------|------------|------------|---------------\n")
  end
  if slv.verbose == 2 || (slv.verbose == 1 && (slv.it == 1 || slv.it%100 == 0))
    @printf("%6d | %7.4e | %7.4e | %11.8e\n", slv.it, slv.gamma, slv.normfpr, slv.cost)
  end
end

function print_status(slv::ForwardBackwardSolver, verbose::Int)
  if verbose == 2 || (verbose == 1 && (slv.it == 1 || slv.it%100 == 0))
    @printf("%6d | %7.4e | %7.4e | %11.8e\n", slv.it, slv.gamma, slv.normfpr, slv.cost)
  end
end

function halt(tol::Float64, gamma::Float64, fpr0::Float64, fpr::Float64, fun_prev::Float64, fun::Float64)
	conv_fpr = fpr <= (1+fpr0)*tol
	conv_fun = abs(fun-fun_prev) <= (1+abs(fun))*tol
	return conv_fpr && conv_fun
end

function Base.show(io::IO, slv::ForwardBackwardSolver)
  println(io, slv.name)
  println(io, "iterations : $(slv.it) / $(slv.maxit)")
  println(io, "fpr        : $(slv.normfpr)")
  println(io, "cost       : $(slv.cost)")
  println(io, "γ          : $(slv.gamma)")
  print(  io, "time       : $(slv.time)")
end
