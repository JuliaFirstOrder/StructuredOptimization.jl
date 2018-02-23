immutable Expression{N} <: AbstractExpression
	x::NTuple{N,Variable}
	L::AbstractOperator
	d::Union{Number, AbstractArray}
	function Expression{N}(x::NTuple{N,Variable},L,d) where {N}

		# checks on L
		ndoms(L,1) > 1 && throw(ArgumentError(
	"cannot create expression with LinearOperator with `ndoms(L,1) > 1`"))

		#checks on x
		szL = size(L,2)
		szx = size.(x)
		check_sz = length(szx) == 1 ? szx[1] != szL : szx != szL
		check_sz && throw(ArgumentError(
	"Size of the operator domain $(size(L, 2)) must match size of the variable $(size.(x))"))

		dmL = domainType(L)
		dmx = eltype.(x)
		check_dm = length(dmx) == 1 ? dmx[1] != dmL : dmx != dmL
		check_dm && throw(ArgumentError(
	"Type of the operator domain $(domainType(L)) must match type of the variable $(eltype.(x))"))

		#checks on d
		if typeof(d) <: AbstractArray
			size(L,1) != size(d) && throw(ArgumentError(
	"cannot sum Array of dimensions $(size(d)) with expression with codomain size $(size(L,1))"))
			codomainType(L) != eltype(d) && throw(ArgumentError(
	"cannot sum $(typeof(d)) with expression with codomain $(codomainType(L))"))
		else
			if typeof(d) <: Complex && codomainType(L) <: Real
				throw(ArgumentError(
	"cannot sum $(typeof(d)) to an expression with codomain $(codomainType(L))"))
			end
		end
		new{N}(x,L,d)
	end
end

include("utils.jl")
include("multiplication.jl")
include("addition.jl")
include("abstractOperator_bind.jl")
