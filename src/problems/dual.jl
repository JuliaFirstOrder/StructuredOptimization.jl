
immutable Dual{Op<:AffineOperator}
	A::Op
	g::ProximableFunction
	y::AbstractArray #dual variables are stored here
end

function problem{T<:NonSmoothTerm}(fi::T, smooth::Array{OptTerm,1}, nonsmooth::Array{OptTerm,1})

	x = fi.A.x #extract variables
	bs = Array{AbstractArray,1}()
	Ainv = Array{LinearOp,1}()
	bi = []
			
	for s in smooth
		if isDiagonal(s.A) == false error("linear operator must invertible") end
		if typeof(s.A) <: Affine 
			push!(bs,s.A.b) 
			push!(Ainv,inv(s.A.A))
		else
			push!(Ainv,inv(s.A))
		end
	end
	if typeof(fi.A) <: Affine bi = copy(fi.A.b) end

	#create dual variables
	y = fi.A*x.x
	#create dual operator
	At = Ainv[1]'*fi.A'
	isempty(bs) ?  nothing : At = At-bs[1]

	g = Precompose(Conjugate(get_prox(fi)),-1,0)
	#is Precompose really needed?
	isempty(bi) ? nothing : g = Tilt(g,bi,0)
	 
	Dual(At, g, y)
	
end

function solve(P::Dual, args...) 
	y, slv = solve!(P.y, P.A, P.g, args...)
	return P.A*y, slv
end


