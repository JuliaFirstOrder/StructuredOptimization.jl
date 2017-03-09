
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

<={K <: AffineOperator, S1 <: Union{Real, AbstractArray}}(A::K,  ub::S1) = IndBox(A,-Inf, ub ) 
<={K <: AffineOperator, S1 <: Union{Real, AbstractArray}}(lb::S1,  A::K) = IndBox(A,  lb, Inf) 


in(x::OptVar, args...) = in(eye(x), args...)
function in{K <: AffineOperator, S1 <: Union{Real, AbstractArray}}(A::K, bnd::AbstractArray{S1,1}) 
	if length(bnd) != 2 error("should provide 2 bounds!") end
	IndBox(A,bnd[1],bnd[2])
end

get_prox(T::IndBox) = ProximalOperators.IndBox(T.lb,T.ub)

fun_name(T::IndBox) = " ⋅ ∈ [l, u] "
fun_par{K,S1<:Real,S2<:Real}(T::IndBox{K,S1,S2}) = " l = $(round(T.lb,3)), u = $(round(T.ub,3)) "
fun_par{K,S1<:Real,S2<:AbstractArray}(T::IndBox{K,S1,S2}) = 
" l = $(round(T.lb,3)), u = $(typeof(T.ub)) "
fun_par{K,S1<:AbstractArray,S2<:Real}(T::IndBox{K,S1,S2}) = 
" l = $(typeof(T.lb)), u = $(round(T.ub,3)) "
fun_par{K,S1<:AbstractArray,S2<:AbstractArray}(T::IndBox{K,S1,S2}) = 
" l = $(typeof(T.lb)), u = $(typeof(T.ub)) "
