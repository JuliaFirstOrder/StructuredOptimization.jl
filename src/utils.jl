module Utils

export print_status

function print_status(it, gamma, fpr, verbose)
  if verbose > 0 && it == 1
    @printf("%6s | %10s | %10s\n", "it", "gamma", "fpr")
    @printf("-------|------------|-----------\n")
  end
  if verbose == 2 || (verbose == 1 && (it == 1 || it%100 == 0))
    @printf("%6d | %7.4e | %7.4e\n", it, gamma, fpr)
  end
end

end
