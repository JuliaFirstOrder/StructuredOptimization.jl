immutable IndBox{S1 <: Union{Real, AbstractArray}, 
		 S2 <: Union{Real, AbstractArray}}   <: NonSmoothFunction
	lb::S1
	ub::S2
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
