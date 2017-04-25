
immutable MoreauEnvelope <: SmoothFunction
	lambda::Real
	p::ProximableFunction
	buf::AbstractVector{Nullable{AbstractArray}}
	# dirty trick to use in place prox! when evaluating the function
	# not sure about that!

	MoreauEnvelope(lambda,p) = new(lambda, p, [ Nullable{AbstractArray}() ])
end

function (s::MoreauEnvelope)(x::AbstractArray)
	if isnull(s.buf[1])
		s.buf[1] = Nullable{AbstractArray}(similar(x))
	end
	prox!(get(s.buf[1]), s.p, x, s.lambda)
	return s.p(get(s.buf[1]))+1/(2*s.lambda)*vecnorm(get(s.buf[1])-x)^2
end

function gradient!(grad::AbstractArray, s::MoreauEnvelope, x::AbstractArray)  
	prox!(grad, s.p, x, s.lambda)
	fx = s.p(grad)
	grad .=  x - grad 
	fx = fx+1/(2*s.lambda)*deepvecnorm(grad)^2
	grad .*= 1/s.lambda
	return fx
end

function gradstep!(x1::AbstractArray, s::MoreauEnvelope, x0::AbstractArray, gamma::Real)  
	gradient!(x1,s,x0)
	x1 .= x0 .- gamma.*x1
	fx1 = s(x1)
	return fx1
end

fun_name(f::MoreauEnvelope,i::Int64) = 
"f$(i)(prox{λ$(i),f$(i)}(A$(i)x))+ 1/2 ‖x - prox{λ$(i),f$(i)}(A$(i)x)‖²"

fun_par( f::MoreauEnvelope,i::Int64)  = "λ$i = $(round(f.lambda,3))"
