export smooth

immutable SmoothedFunction <: SmoothFunction
	gamma::Real
	p::ProximableFunction
	f::NonSmoothFunction
end

function (s::SmoothedFunction)(x::AbstractArray)
	return s.p(x)+0.5*s.gamma*vecnorm(x)^2
end

function gradient!(grad::AbstractArray, s::SmoothedFunction, x::AbstractArray)  
	prox!(s.p,x,grad,s.gamma)
	grad .= 1/s.gamma.*( x - grad )
end

function smooth(cf::CostFunction, gamma::Real=1.)
	f = Vector{ExtendedRealValuedFunction}(0)
	for fi in terms(cf)
		if typeof(fi)<:SmoothFunction
			push!(f,fi)
		elseif typeof(fi)<:NonSmoothFunction
			p = get_prox(fi)
			fs = SmoothedFunction(gamma,p,fi) 
			push!(f,fs)
		end
	end
	CostFunction(variable(cf),f,affine(cf))
end

fun_name(f::SmoothedFunction,i::Int64) =  fun_name(f.f,i)*"+ γ$i/2 ‖A$(i)x‖²"
fun_par( f::SmoothedFunction,i::Int64)  = fun_par(f.f,i)*", γ$i = $(round(f.gamma,3))"
