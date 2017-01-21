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

function halt(slv::ForwardBackwardSolver, normfpr0::Float64, cost_prev::Float64)
	conv_fpr = slv.normfpr <= (1+normfpr0)*slv.tol
	conv_fun = abs(slv.cost-cost_prev) <= (1+abs(slv.cost))*slv.tol
	return conv_fpr && conv_fun
end

function halt(slv::ZeroFPR, normfpr0::Float64, FBEx::Float64, FBEprev::Float64,)
	conv_fpr = slv.normfpr <= (1+normfpr0)*slv.tol
	conv_fun = abs(FBEx-FBEprev) <= (1+abs(FBEx))*slv.tol
	return conv_fpr && conv_fun
end

# compute upper bound for Lipschitz constant using finite differences
function get_gamma0{T<:Union{Complex{Float64},Float64}}(L::Function, Ladj::Function, x::AbstractArray{T}, gradx::AbstractArray, b::AbstractArray)
	resy  = L(x+sqrt(eps())) - b
	grady = Ladj(resy)
	return vecnorm( sqrt(eps())*ones(x))/vecnorm(gradx-grady)
end

function get_gamma0{T<:AbstractArray}(L::Function, Ladj::Function, x::AbstractArray{T}, gradx::AbstractArray, b::AbstractArray)
	resy  = L(x+sqrt(eps())) - b
	grady = Ladj(resy)
	z = similar(x)
	[z[i] = ones(x[i]) for i in eachindex(z)]
	return deepvecnorm( sqrt(eps())*z)/deepvecnorm(gradx-grady)
end

deepvecnorm{T<:Union{Complex{Float64},Float64}}(x::AbstractArray{T}) = vecnorm(x)

function deepvecnorm{T<:AbstractArray}(x::AbstractArray{T})
	out = 0.
	for a in x
		out += vecnorm(a)^2
	end
	return sqrt(out)
end


function Base.show(io::IO, slv::ForwardBackwardSolver)
  println(io, slv.name)
  println(io, "iterations : $(slv.it) / $(slv.maxit)")
  println(io, "fpr        : $(slv.normfpr)")
  println(io, "cost       : $(slv.cost)")
  println(io, "γ          : $(slv.gamma)")
  println(io, "time       : $(slv.time)")
end
