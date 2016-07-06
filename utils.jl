function print_status()
  @printf("%6s | %10s | %10s\n", "it", "gamma", "fpr")
  @printf("-------|------------|-----------\n")
end

function print_status(it, gamma, fpr)
  @printf("%6d | %7.4e | %7.4e\n", it, gamma, fpr)
end
