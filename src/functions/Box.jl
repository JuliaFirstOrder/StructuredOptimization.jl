
import Base: <=, in 

type IndBox{K <: AffineOperator, 
	    S1 <: Union{Real, AbstractArray}, 
	    S2 <: Union{Real, AbstractArray}}   <: NonSmoothTerm
	A::K
	lb::S1
	ub::S2
end

<=(x::OptVar, args...) = <=(eye(x), args...)
<={S1 <: Union{Real, AbstractArray}}(lb::S1, x::OptVar) = <=(lb, eye(x))

<={K <: AffineOperator, S1 <: Union{Real, AbstractArray}}(A::K,  ub::S1) = IndBox(A,-Inf, ub) 
<={K <: AffineOperator, S1 <: Union{Real, AbstractArray}}(lb::S1,  A::K) = IndBox(A,  lb,Inf) 

in(x::OptVar, args...) = in(eye(x), args...)
function in{K <: AffineOperator, S1 <: Union{Real, AbstractArray}}(A::K, bnd::AbstractArray{S1,1}) 
	if length(bnd) != 2 error("should provide 2 bounds!") end
	IndBox(A,bnd[1],bnd[2])
end

function get_prox(T::IndBox)
	return ProximalOperators.IndBox(T.lb,T.ub)
end
