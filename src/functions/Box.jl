
import Base: <=, in 

type IndBox{S1 <: Union{Real, AbstractArray}, 
	    S2 <: Union{Real, AbstractArray}}   <: NonSmoothFunction
	lb::S1
	ub::S2
end

<=(x::OptVar, args...) = <=(eye(x), args...)
<={S1 <: Union{Real, AbstractArray}}(lb::S1, x::OptVar) = <=(lb, eye(x))

<={S1 <: Union{Real, AbstractArray}}(A::AffineOperator,  ub::S1) =  
CostFunction(variable(A), [IndBox(-Inf, ub )], [A])
<={S1 <: Union{Real, AbstractArray}}(lb::S1,  A::AffineOperator) =  
CostFunction(variable(A), [IndBox(lb, Inf )], [A])


in(x::OptVar, args...) = in(eye(x), args...)
function in{S1 <: Union{Real, AbstractArray}}(A::AffineOperator, bnd::AbstractArray{S1,1}) 
	if length(bnd) != 2 error("should provide 2 bounds!") end
	CostFunction(variable(A), [IndBox(bnd[1],bnd[2])], [A])
end

get_prox(T::IndBox) = ProximalOperators.IndBox(T.lb,T.ub)

fun_name(T::IndBox,i::Int64) = "Ind{A$(i)x ∈ [l$(i), u$(i)]}(x)"
fun_par{S1<:Real,S2<:Real}(T::IndBox{S1,S2},i::Int64) = 
" l$(i) = $(round(T.lb,3)), u$(i) = $(round(T.ub,3)) "
fun_par{S1<:Real,S2<:AbstractArray}(T::IndBox{S1,S2},i::Int64) = 
" l$(i) = $(round(T.lb,3)), u$(i) = $(typeof(T.ub)) "
fun_par{S1<:AbstractArray,S2<:Real}(T::IndBox{S1,S2},i::Int64) = 
" l$(i) = $(typeof(T.lb)), u$(i) = $(round(T.ub,3)) "
fun_par{S1<:AbstractArray,S2<:AbstractArray}(T::IndBox{S1,S2},i::Int64) = 
" l$(i) = $(typeof(T.lb)), u$(i) = $(typeof(T.ub)) "
