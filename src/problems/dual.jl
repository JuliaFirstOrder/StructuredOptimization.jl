
immutable Dual
	x::Vector{Variable} #primal variables
	cf::CostFunction
	p::ProximableFunction
	Ainv::Vector{LinearOperator}
end

function Dual(smooth::CostFunction, nonsmooth::CostFunction)

	length(terms(smooth)) != length(variable(smooth)) && error("not enough terms smooth for dual") 
	!(isLeastSquares(smooth)) && error("terms in smooth must all be least squares")

	y = Variable(affine(nonsmooth)[1](~variable(nonsmooth)))     #create dual variables
	B = [operator(affine(nonsmooth)[1])]     #extract operator form nonsmooth
	typeof(B[1]) <: HCAT ? B = B[1].A : nothing

	Ainv = Array{LinearOperator,1}(length(terms(smooth)))

	for i in eachindex(variable(smooth))
		A = operator(affine(smooth)[i])
		typeof(A)<: HCAT ? A = A.A[i] : nothing
		isInvertable(A) == true ? nothing : error("operator not easily invertable")
		Ainv[i] = inv(A)
	end
		
	cf = CostFunction() #new dual cost function
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
	Precompose(Conjugate(get_prox(g)),-1,0)
end

function solve(P::Dual, slv::Solver = default_slv()) 
	slv = solve(P.cf, P.p, slv)
	resy, = evaluate(P.cf,~variable(P.cf)) 
	~P.x[1] .= -(P.Ainv[1]'*resy[1])
	for i = 2:length(resy)
		~P.x[i] .= -(P.Ainv[i]'*resy[i])
	end
	return slv
end


function Base.show(io::IO, P::Dual)
	println("Dual Problem")
	println()
	println("Smooth Cost Function:")
	println()
	show(P.cf)
	println()
	println("Proximable operators:")
	println()
	show(P.p)
	println()
end

