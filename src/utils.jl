function print_status(it, gamma, fpr, cost, verbose)
  if verbose > 0 && it == 1
    @printf("%6s | %10s | %10s | %14s\n", "it", "gamma", "fpr", "cost")
    @printf("-------|------------|------------|---------------\n")
  end
  if verbose == 2 || (verbose == 1 && (it == 1 || it%100 == 0))
    @printf("%6d | %7.4e | %7.4e | %11.8e\n", it, gamma, fpr, cost)
  end
end
