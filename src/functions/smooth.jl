export smooth

immutable MoreauEnvelope <: SmoothFunction
	gamma::Real
	p::ProximableFunction
	f::NonSmoothFunction
end

function (s::MoreauEnvelope)(x::AbstractArray)
	return s.p(x)+0.5*s.gamma*vecnorm(x)^2
end

function gradient!(grad::AbstractArray, s::MoreauEnvelope, x::AbstractArray)  
	prox!(s.p,x,grad,s.gamma)
	grad .= 1/s.gamma.*( x - grad )
end

function smooth(cf::CostFunction, gamma0::Real=1.)
	f = Vector{ExtendedRealValuedFunction}(0)
	for fi in terms(cf)
		if typeof(fi)<:SmoothFunction
			push!(f,fi)
		elseif typeof(fi)<:NonSmoothFunction
			p = get_prox(fi)
			fs = MoreauEnvelope(gamma0*lambda(fi),p,fi) 
			push!(f,fs)
		end
	end
	CostFunction(variable(cf),f,affine(cf))
end

fun_name(f::MoreauEnvelope,i::Int64) =  fun_name(f.f,i)*"+ γ$i/2 ‖A$(i)x‖²"
fun_par( f::MoreauEnvelope,i::Int64)  = fun_par(f.f,i)*", γ$i = $(round(f.gamma,3))"
