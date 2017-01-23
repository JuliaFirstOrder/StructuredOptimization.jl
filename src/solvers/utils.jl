# Solver trace

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

# Stopping criteria

function halt(slv::ForwardBackwardSolver, normfpr0::Float64, cost_prev::Float64)
	conv_fpr = slv.normfpr <= (1+normfpr0)*slv.tol
	conv_fun = abs(slv.cost-cost_prev) <= (1+abs(slv.cost))*slv.tol
	return conv_fpr && conv_fun
end

function halt(slv::ZeroFPR, normfpr0::Float64, FBEx::Float64, FBEprev::Float64)
	conv_fpr = slv.normfpr <= (1+normfpr0)*slv.tol
	conv_fun = abs(FBEx-FBEprev) <= (1+abs(FBEx))*slv.tol
	return conv_fpr && conv_fun
end

# Generalized length, dot product and norm, for nested Array objects

function deeplength(x::AbstractArray)
  len = 0
	for k in eachindex(x)
		len += deeplength(x[k])
	end
	return len
end

deeplength{T <: Number}(x::AbstractArray{T}) = length(x)

function deepvecdot(x::AbstractArray, y::AbstractArray)
	out = 0.0
	for k in eachindex(x)
		out += deepvecdot(x[k], y[k])
	end
	return out
end

deepvecdot{T <: Number}(x::AbstractArray{T}, y::AbstractArray{T}) = vecdot(x, y)

deepvecnorm(x::AbstractArray) = sqrt(deepvecdot(x, x))

deepvecnorm{T <: Number}(x::AbstractArray{T}) = vecnorm(x)

# To print out solver objects

function Base.show(io::IO, slv::ForwardBackwardSolver)
  println(io, slv.name)
  println(io, "iterations : $(slv.it) / $(slv.maxit)")
  println(io, "fpr        : $(slv.normfpr)")
  println(io, "cost       : $(slv.cost)")
  println(io, "γ          : $(slv.gamma)")
  println(io, "time       : $(slv.time)")
end
