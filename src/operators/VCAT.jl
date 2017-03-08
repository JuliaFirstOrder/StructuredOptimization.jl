import Base: vcat

type VCAT{D3} <: LinearOperator{D3}
	x::Array{OptVar}
	A::AbstractArray{LinearOperator}
	mid::AbstractArray
	sign::Array{Bool,1}
end

function size(A::VCAT)
	dim1, dim2 = size(A.A[1],2), tuple(size.(A.A)[2]...)  
end

fun_name(S::VCAT) = "Vertically Concatenated Operators"

function vcat{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3}, sign::Array{Bool} )
	if size(A,1) != size(B,1) DimensionMismatch("operators must share codomain!") end
	mid = Array{D3}(size(B,1))
	VCAT{D3}([A.x], [A,B], mid, sign)
end

vcat{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3} ) = vcat(A,B,[true; true])
vcat(A::LinearOperator, B::OptVar  ) = vcat(A     , eye(B))
vcat(A::OptVar  , B::LinearOperator) = vcat(eye(A),     B )
vcat(A::OptVar  , B::OptVar  ) = vcat(eye(A), eye(B))

function vcat(A::Vararg{Union{LinearOperator,OptVar}})
	V = vcat(A[1],A[2])
	for i = 3:length(A)
		typeof(A[i]) <: OptVar ? push!(V.A,eye(A[i])) : push!(V.A,A[i])
		push!(V.sign,true)
	end
	return V
end

transpose{D3}(A::VCAT{D3}) = HCAT{D3}(A.x, A.A.', A.mid, A.sign)

function *{D2}(A::VCAT{D2},b::AbstractArray) 
	y = Array{AbstractArray,1}(length(A.A))
	for i = 1:length(A.A)
		y[i] = create_out(A.A[i]) 
	end
	A_mul_B!(y,A,b)
	return y
end

function A_mul_B!{T1<:AbstractArray}(y::Array{T1,1}, S::VCAT, b::AbstractArray) 
	for i = 1:length(S.A)
		S.sign[i] ? A_mul_B!(y[i],S.A[i],b) : A_mul_B!(y[i],S.A[i],-b)
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

