type CostFunction
	x::Array{OptVar,1}
	Terms::Array{OptTerm, 1}
end

+(h::OptTerm, g::OptTerm) = CostFunction(addVar(variable(h),variable(g)),[h,g])

function +(cf::CostFunction, g::OptTerm) 
	push!(cf.Terms,g)
	addVar!(cf.x,variable(g))
	return cf
end

function addVar{T1,T2}(x::OptVar{T1},y::OptVar{T2})
	x == y ? [x] :[x,y]
end

function addVar!{T}(x::Vector,y::OptVar{T})
	any(x.==y) ? x : push!(x,y) 
end
addVar!{T}(y::OptVar{T}, x::Vector) = addVar!(x,y) 

function addVar{T}(x::Vector,y::OptVar{T})
	any(x.==y) ? x : [x...,y] 
end
addVar{T}(y::OptVar{T}, x::Vector) = addVar(x,y) 

function addVar!(x::Vector,y::Vector)
	for xi in x
		addVar!(xi,z)
	end
end

function addVar(x::Vector,y::Vector)
	z = copy(y)
	for xi in x
		z = addVar(xi,z)
	end
	return z
end
