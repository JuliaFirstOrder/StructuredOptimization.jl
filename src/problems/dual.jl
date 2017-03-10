
immutable Dual
	At::AffineOperator
	Ainv::LinearOperator
	p::ProximableFunction
	y::AbstractArray #dual variables are stored here
end

function problem{T<:NonSmoothTerm}(g::T, smooth::Array{OptTerm,1})

	x = variable(g)        #extract variables
	# TODO check that length(smooth) == number of blocks of variables
	y = g.A*optArray(g.A)     #create dual variables

	b2 = 0.*copy(y)
	if typeof(g.A) <: Affine b2 .-= g.A.b end

	Ainv = Array{LinearOperator,1}(length(smooth))
	b1 = 0.*deepcopy(optArray(g.A))
	for s in smooth
		isInvertable(s.A) == true ? nothing : error("operator not easily invertable")
		if typeof(s.A) <: Affine
			get_b1!(b1,x,s)
			get_Ainv!(Ainv,x,s)
		else
			get_Ainv!(Ainv,x,s)
		end
	end

	Ainv = hcat(Ainv...)
	At = Ainv'.*g.A'
	At = At+b1

	p = get_cjprox(g,b2) 
	 
	Dual(At, Ainv, p, y)
	
end

function get_b1!(b1::AbstractArray, x::OptVar, s::LeastSquares)
	b1 .= s.A.b
end

function get_b1!{T<:AbstractArray}(b1::Array{T,1}, x::Array{OptVar,1}, s::LeastSquares)
	for i in eachindex(x)
		if s.A.x == x[i]
			b1[i] .= s.A.b
		end
	end
end

function get_Ainv!(Ainv::Array{LinearOperator,1}, x::OptVar, s::LeastSquares)
	Ainv[1] = inv(s.A)
end

function get_Ainv!(Ainv::Array{LinearOperator,1}, x::Array{OptVar,1}, s::LeastSquares)
	for i in eachindex(x)
		if s.A.x == x[i]
			Ainv[i] = inv(s.A)
		end
	end
end

function get_cjprox{T<:NonSmoothTerm}(g::T, b2::AbstractArray)
	#is precompose really needed?
	Tilt(Precompose(Conjugate(get_prox(g)),-1,0), b2, 0.0)
end

function solve{S<:Solver}(P::Dual, slv::S) 
	slv2 = copy(slv)
	y, slv2 = solve!(P.y, P.At, P.p, slv2)
	return -(P.Ainv.*(P.At*y)), slv2
end


