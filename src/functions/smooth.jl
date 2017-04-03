export smooth

immutable MoreauEnvelope <: SmoothFunction
	gamma::Real
	p::ProximableFunction
	f::NonSmoothFunction
end

function (s::MoreauEnvelope)(x::AbstractArray)
	p, = prox(s.p, x, s.gamma)
	return s.p(p)+1/(2*s.gamma)*vecnorm(p-x)^2
end

function gradient!(grad::AbstractArray, s::MoreauEnvelope, x::AbstractArray)  
	prox!(s.p, x, grad, s.gamma)
	fgamma = s.p(grad)
	grad .=  x - grad 
	fx = fgamma+1/(2*s.gamma)*vecnorm(grad)^2
	grad .*= 1/s.gamma
	return fx
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

fun_name(f::MoreauEnvelope,i::Int64) =  "f$(i)(prox{γ$(i),f$(i)}(A$(i)x))+ 1/2 ‖x - prox{γ$(i),f$(i)}(A$(i)x)‖²"
fun_par( f::MoreauEnvelope,i::Int64)  = fun_par(f.f,i)*", γ$i = $(round(f.gamma,3))"
