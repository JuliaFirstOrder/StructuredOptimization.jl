import Base: hcat

type HCAT{D3} <: LinearOperator{D3}
	A::AbstractArray{LinearOperator}
	mid::AbstractArray
	sign::Array{Bool,1}
end

variable(A::HCAT) = variable.(A.A)

function size(A::HCAT)
	dim1, dim2 = tuple(size.(A.A)[2]...), size(A.A[1],2)   
end

fun_name(S::HCAT) = "Horizontally Concatenated Operators"

hcat(A::LinearOperator) = A

function hcat{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3}, sign::Array{Bool} )
	if size(A,2) != size(B,2) DimensionMismatch("operators must go to same space!") end
	mid = Array{D3}(size(B,2))
	HCAT{D3}([A,B], mid, sign)
end

hcat{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3} ) = hcat(A,B,[true; true])
hcat(A::LinearOperator, B::OptVar  ) = vcat(A     , eye(B))
hcat(A::OptVar  , B::LinearOperator) = vcat(eye(A),     B )
hcat(A::OptVar  , B::OptVar  ) = vcat(eye(A), eye(B))

function hcat(A::Vararg{Union{LinearOperator,OptVar}})
	H = hcat(A[1],A[2])
	for i = 3:length(A)
		typeof(A[i]) <: OptVar ? push!(H.A,eye(A[i])) : push!(H.A,A[i])
		push!(H.sign,true)
	end
	return H
end

transpose{D3}(A::HCAT{D3}) = VCAT{D3}(A.A.', A.mid, A.sign)

function *{D3,T1<:AbstractArray}(A::HCAT{D3},b::Array{T1,1}) 
	y = Array{D3}(size(A,2))
	A_mul_B!(y,A,b)
	return y
end

function .*{D2,T1<:AbstractArray}(A::HCAT{D2},b::Array{T1,1}) 
	y = Array{AbstractArray,1}(length(A.A))
	for i = 1:length(A.A)
		y[i] = create_out(A.A[i]) 
		A_mul_B!(y[i],A.A[i],b[i])
	end
	return y
end

.*{D3}(A::HCAT{D3},B::HCAT{D3})  = HCAT{D3}(A.A.*B.A, A.mid, A.sign.*B.sign)

#forward
function A_mul_B!{T1<:AbstractArray}(y::AbstractArray,S::HCAT,b::Array{T1,1}) 
	S.sign[1] ? A_mul_B!(y,S.A[1],b[1]) : A_mul_B!(y,S.A[1],-b[1])
	for i = 2:length(S.A)
		A_mul_B!(S.mid,S.A[i],b[i])
		S.sign[i] ? y .= (+).(y,S.mid) : y .= (-).(y,S.mid)
	end
end

create_out{D1,D2}(A::LinearOperator{D1,D2}) = Array{D2}(size(A,2))

#adjoint
function A_mul_B!{T1<:AbstractArray}(y::Array{T1,1},S::HCAT,b::AbstractArray) 
	for i = 1:length(S.A)
		S.sign[i] ? A_mul_B!(y[i],S.A[i],b) : A_mul_B!(y[i],S.A[i],-b)
	end
end

#printing stuff
function fun_dom{D3<:Real}(A::HCAT{D3}) 
	str = ""
	for a in A.A str *= fun_D1(a,1) end
	str *= "→  ℝ^$(size(A,2))"
end

function fun_dom{D3<:Complex}(A::HCAT{D3}) 
	str = ""
	for a in A.A str *= fun_D1(a,1) end
	str *= "→  ℂ^$(size(A,2))"
end

fun_D1{D1<:Real, D2}(A::LinearOperator{D1,D2},dim::Int64)    =  " ℝ^$(size(A,dim)) "
fun_D1{D1<:Complex, D2}(A::LinearOperator{D1,D2},dim::Int64) =  " ℂ^$(size(A,dim)) "

import Base: copy, sortperm, sort!, sort, issorted

copy{D3}(A::HCAT{D3}) = HCAT{D3}(copy(A.A),A.mid,copy(A.sign))

sortperm(A::HCAT) = sortperm(variable(A), by = object_id)
sortperm(A::LinearOperator) = false

sortperm(A::Affine) = sortperm(A.A) 

function sort!(A::HCAT)
	p = sortperm(A) 
	A.A .= A.A[p]
	A.sign .= A.sign[p]
end

function sort(A::HCAT)
	A2 = copy(A)
	sort!(A2)
	return A2
end

function sort(A::Affine)
	p = sortperm(A)
	if p == false 
		return A
	else #it's HCAT
		A2 = sort(A.A)
		return A2+A.b
	end
end

issorted(A::HCAT) = issorted(sortperm(A))

sort!(A::LinearOperator) = nothing
sort!(A::Affine) = sort!(A.A)
sort(A::LinearOperator)  = A
issorted(A::LinearOperator) = true
issorted(A::Affine) = issorted(A.A)



