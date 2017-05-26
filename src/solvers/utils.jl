# Solver trace

function print_status(slv::ForwardBackwardSolver)
  if slv.verbose > 0 && slv.it == 1
    @printf("%6s | %10s | %10s | %14s\n", "it", "gamma", "fpr", "cost")
    @printf("-------|------------|------------|---------------\n")
  end
  if slv.verbose == 2 || (slv.verbose == 1 && (slv.it == 1 || slv.it%100 == 0))
    @printf("%6d | %7.4e | %7.4e | %11.8e\n", slv.it, slv.gamma, slv.normfpr/slv.gamma, slv.cost)
  end
end

function print_status(slv::ForwardBackwardSolver, verbose::Int)
  if verbose == 2 || (verbose == 1 && (slv.it == 1 || slv.it%100 == 0))
    @printf("%6d | %7.4e | %7.4e | %11.8e\n", slv.it, slv.gamma, slv.normfpr/slv.gamma, slv.cost)
  end
end

# Stopping criteria

function halt_default(slv::ForwardBackwardSolver)
	conv_fpr = slv.normfpr/slv.gamma <= slv.tol
	return conv_fpr
end

# function halt_default(slv::ZeroFPR, normfpr0::Float64, FBEx::Float64, FBEprev::Float64)
# 	conv_fpr = slv.normfpr <= (1+normfpr0)*slv.tol
# 	conv_fun = abs(FBEx-FBEprev) <= (1+abs(FBEx))*slv.tol
# 	return conv_fpr && conv_fun
# end
