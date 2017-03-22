import Base: hcat

immutable HCAT{D3} <: LinearOperator{D3}
	A::Vector{LinearOperator}
	mid::AbstractArray
	sign::Array{Bool,1}
end

  domainType(A::HCAT) = domainType.(A.A)
codomainType{D3}(A::HCAT{D3}) = D3

function size(A::HCAT)
	dim1, dim2 = tuple(size.(A.A)[2]...), size(A.A[1],2)   
end

fun_name(S::HCAT) = "Horizontally Concatenated Operators"

hcat(A::LinearOperator) = A

function hcat{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3}, sign::Array{Bool} )
	if size(A,2) != size(B,2) throw(DimensionMismatch("operators must go to same space!")) end
	mid = Array{D3}(size(B,2))
	HCAT{D3}([A,B], mid, sign)
end

hcat{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3} ) = hcat(A,B,[true; true])

function hcat(A::Vararg{LinearOperator})
	H = Vector{LinearOperator}()
	D3   = codomainType(A[1])
	dim2 = size(A[1],2)
	mid = Array{D3}(dim2)
	for a in A
		if size(a,2) != dim2 || codomainType(a) != D3
			throw(DimensionMismatch("operators must go to same space!"))
		end
		push!(H,a)
	end
	sign = ones(Bool,length(A))
	HCAT{D3}(H,mid,sign)
end

transpose{D3}(A::HCAT{D3}) = VCAT{D3}((A.A.')[:], A.mid, A.sign)

function *{D3,T1<:AbstractArray}(A::HCAT{D3},b::Array{T1,1}) 
	y = Array{D3}(size(A,2))
	A_mul_B!(y,A,b)
	return y
end

function .*{D2,T1<:AbstractArray}(A::HCAT{D2},b::Array{T1,1}) 
	y = Array{AbstractArray,1}(length(A.A))
	for i = 1:length(A.A)
		y[i] = A.A[i]*b[i] 
	end
	return y
end

.*{D3}(A::HCAT{D3},B::HCAT{D3})  = HCAT{D3}(A.A.*B.A, A.mid, A.sign.*B.sign)

function A_mul_B!{T1<:AbstractArray}(y::AbstractArray,S::HCAT,b::Array{T1,1}) 
	S.sign[1] ? A_mul_B!(y,S.A[1],b[1]) : A_mul_B!(y,S.A[1],-b[1])
	for i = 2:length(S.A)
		A_mul_B!(S.mid,S.A[i],b[i])
		S.sign[i] ? y .= (+).(y,S.mid) : y .= (-).(y,S.mid)
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

import Base: copy, sort

copy{D3}(A::HCAT{D3}) = HCAT{D3}(copy(A.A),A.mid,copy(A.sign))

function sort{D3}(A::HCAT{D3},p::Array)
	H = A.A[p]
	s = A.sign[p]
	return HCAT{D3}(H,A.mid,s)
end


extract_operator(A::HCAT          , idx::Int64) = A.A[idx], A.sign[idx]


