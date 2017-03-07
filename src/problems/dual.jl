
immutable Dual
	At::LinearOp
	Ainv::LinearOp
	b::AbstractArray
	g::ProximableFunction
	y::AbstractArray #dual variables are stored here
end

function problem{T<:NonSmoothTerm}(fi::T, smooth::Array{OptTerm,1})

	x = fi.A.x #extract variables
	typeof(smooth[1].A) <: Affine ? Ainv = inv(smooth[1].A.A) : Ainv = inv(smooth[1].A)
	At = Ainv*fi.A'


	y = fi.A*fi.A.x.x  #create dual variables
	b = deepsimilar(y)

	if typeof(fi.A) <: Affine        b -= fi.A.b end
	if typeof(smooth[1].A) <: Affine b += At'*smooth[1].A.b end

	g = Tilt(Precompose(Conjugate(get_prox(fi)),-1,0), b, 0.0)
	 
	if typeof(smooth[1].A) <: Affine
		Dual(At, Ainv, smooth[1].A.b, g, y)
	else
		Dual(At, Ainv, [], g, y)
	end
	
end

function solve{S<:Solver}(P::Dual, slv::S) 
	slv2 = copy(slv)
	y, slv2 = solve!(P.y, P.At, P.g, slv2)
	if isempty(P.b)
		return -(P.Ainv*(P.At*y)), slv2
	else
		return -(P.Ainv*(P.At*y+P.b)), slv2
	end
end


