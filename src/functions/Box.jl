
import Base: <= 

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
#<={S1 <: Union{Real, AbstractArray}}(lb::S1,A::IndBox) = IndBox(A.A,lb, A.ub) 
<={S1 <: Union{Real, AbstractArray}}(A::IndBox,ub::S1) = IndBox(A.A,A.lb,ub) 

