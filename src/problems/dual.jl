
immutable Dual
	x::Vector{Variable} #primal variables
	cf::Term
	p::ProximableFunction
	Ainv::Vector{AbstractOperator}
end

function Dual(smooth::Term, nonsmooth::Term)

	length(terms(smooth)) != length(variable(smooth)) && error("not enough terms smooth for dual")
	!(isLeastSquares(smooth)) && error("terms in smooth must all be least squares")

	# create dual variables
	y = Variable(affine(nonsmooth)[1](~variable(nonsmooth)))
	# extract operator form nonsmooth
	B = [operator(affine(nonsmooth)[1])]
	typeof(B[1]) <: HCAT ? B = B[1].A : nothing

	Ainv = Array{AbstractOperator,1}(length(terms(smooth)))

	for i in eachindex(variable(smooth))
		A = operator(affine(smooth)[i])
		typeof(A)<: HCAT ? A = A.A[i] : nothing
		is_invertible(A) == true ? nothing : error("operator not easily invertable")
		Ainv[i] = inv(A)
	end

	cf = Term() #new dual cost function
	for i in eachindex(variable(smooth))
		lambda = (terms(smooth)[i].lambda)
		if lambda == 1.0
			cf += ls((Ainv[i]'*B[i]')*y + tilt(affine(smooth)[i]) )
		else
			cf += lambda*ls((1/lambda*Ainv[i]'*B[i]')*y + tilt(affine(smooth)[i]) )
		end
	end

	p = get_cjprox(terms(nonsmooth)[1],-tilt(affine(nonsmooth)[1]))

	Dual(variable(smooth),cf, p, Ainv)

end

function get_cjprox{T<:NonSmoothFunction}(g::T, b2::AbstractArray)
	Tilt(get_cjprox(g,0.0), b2, 0.0)
end

#no tilt
function get_cjprox{T<:NonSmoothFunction}(g::T, b2::Float64)
	#is precompose really needed?
	PrecomposeDiagonal(Conjugate(get_prox(g)),-1,0)
end

function solve(P::Dual, slv::Solver = default_slv())
	slv = solve(P.cf, P.p, slv)
	resy, = residual(P.cf,~variable(P.cf))
	~P.x[1] .= -(P.Ainv[1]'*resy[1])
	for i = 2:length(resy)
		~P.x[i] .= -(P.Ainv[i]'*resy[i])
	end
	return slv
end


function Base.show(io::IO, P::Dual)
	println(io, "Dual Problem")
	println(io)
	println(io, "Smooth Cost Function:")
	show(io, P.cf)
	println(io)
	println(io, "Proximable operators:")
	println(io)
	show(io, P.p)
end
