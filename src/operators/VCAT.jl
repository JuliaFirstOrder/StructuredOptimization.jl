import Base: vcat

type VCAT{D3} <: LinearOperator{D3}
	A::AbstractArray{LinearOperator}
	mid::AbstractArray
end
sign(S::VCAT) = true
-{D3}(A::VCAT{D3}) = VCAT{D3}((-).(A.A), A.mid) 

  domainType{D3}(A::VCAT{D3}) = D3
codomainType(A::VCAT) = codomainType.(A.A)

function size(A::VCAT)
	dim1, dim2 = size(A.A[1],2), tuple(size.(A.A)[2]...)  
end

fun_name(S::VCAT) = "Vertically Concatenated Operators"

vcat(A::LinearOperator) = A

function vcat{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3} )
	if size(A,1) != size(B,1) throw(DimensionMismatch("operators must share codomain!")) end
	mid = Array{D3}(size(B,1))
	VCAT{D3}([A,B], mid)
end

function vcat(A::Vararg{LinearOperator})
	H = Vector{LinearOperator}()
	D3   = domainType(A[1])
	dim1 = size(A[1],1)
	mid = Array{D3}(dim1)
	for a in A
		if size(a,1) != dim1 || domainType(a) != D3
			throw(DimensionMismatch("operators must share codomain!"))
		end
		push!(H,a)
	end
	VCAT{D3}(H,mid)
end

transpose{D3}(A::VCAT{D3}) = HCAT{D3}((A.A.')[:], A.mid)

function *{D2}(A::VCAT{D2},b::AbstractArray) 
	y = Array{AbstractArray,1}(length(A.A))
	for i = 1:length(A.A)
		y[i] = Array{D2}(size(A.A[i],2)) 
	end
	A_mul_B!(y,A,b)
	return y
end

.*{D3}(A::VCAT{D3},B::VCAT{D3})  = VCAT{D3}(A.A.*B.A, A.mid)

function A_mul_B!{T1<:AbstractArray}(y::Array{T1,1}, S::VCAT, b::AbstractArray) 
	for i = 1:length(S.A)
		A_mul_B!(y[i],S.A[i],b)
	end
end

#printing stuff
function fun_dom{D3<:Real}(A::VCAT{D3}) 
	str = " ℝ^$(size(A,1)) → "
	for a in A.A str *= fun_D1(a,2) end
	return str
end

function fun_dom{D3<:Complex}(A::VCAT{D3}) 
	str = " ℂ^$(size(A,1)) → "
	for a in A.A str *= fun_D1(a,2) end
	return str
end

