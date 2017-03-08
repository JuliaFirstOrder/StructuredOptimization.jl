
immutable Dual
	At::LinearOperator
	Ainv::LinearOperator
	b::AbstractArray
	g::ProximableFunction
	y::AbstractArray #dual variables are stored here
end

function problem{T<:NonSmoothTerm}(g::T, smooth::Array{OptTerm,1})

	x = g.A.x #extract variables
	isInvertable(smooth[1].A) == true ? nothing : error("operator not easily invertable")
	typeof(smooth[1].A) <: Affine ? Ainv = inv(smooth[1].A.A) : Ainv = inv(smooth[1].A)
	At = g.A'
	At = Ainv'*At

	y = g.A*g.A.x.x  #create dual variables
	b = zeros(y)

	if typeof(g.A) <: Affine         b -= g.A.b end
	if typeof(smooth[1].A) <: Affine b += At'*smooth[1].A.b end

	g = Tilt(Precompose(Conjugate(get_prox(g)),-1,0), b, 0.0)
	 
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


