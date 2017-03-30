import Base: zeros
export Zeros

immutable Zeros{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	dim::Tuple

	Zeros(dim) = new(true,dim)
	Zeros(sign,dim) = new(sign,dim)
end
size(A::Zeros) = (A.dim[1],A.dim[2])
-{D1,D2}(A::Zeros{D1,D2}) = Zeros{D1,D2}(false == sign(A),A.dim) 

Zeros(dim1::Tuple, dim2::Tuple) = Zeros{Float64,Float64}((dim1,dim2))
Zeros(T::Type, dim1::Tuple, dim2::Tuple) = Zeros{T,T}(dim1, dim2)
Zeros(D1::Type, D2::Type, dim1::Tuple, dim2::Tuple) = Zeros{D1,D2}((dim1,dim2))

zeros{D1,D2}(A::LinearOperator{D1,D2}) = Zeros{D1,D2}(size(A,1),size(A,2))

function zeros{D3}(A::VCAT{D3}) 
	dim2 = size(A,2)
	V = Vector{LinearOperator}(length(dim2))
	for i in eachindex(V)
		V[i] = zeros(A.A[i])
	end
	vcat(V...)
end

function zeros{D3}(A::HCAT{D3}) 
	dim1 = size(A,1)
	V = Vector{LinearOperator}(length(dim1))
	for i in eachindex(V)
		V[i] = zeros(A.A[i])
	end
	hcat(V...)
end

zeros{D1}(x::Variable{D1}) = Affine([x], Zeros{D1,D1}((size(x),size(x))), 
				 Zeros{D1,D1}((size(x),size(x))),
				 Nullable{AbstractArray}() )

zeros{D1}(x::Variable{D1}, dim::Vararg{Int64}) = Affine([x], Zeros{D1,D1}((size(x),dim)), 
						     Zeros{D1,D1}((dim,size(x))),
						     Nullable{AbstractArray}() )

transpose{D1}(A::Zeros{D1,D1} ) = Zeros{D1,D1}(sign(A),(A.dim[2],A.dim[1]))

function uA_mul_B!(y::AbstractArray,A::Zeros,b::AbstractArray)
	y .= 0
end

fun_name(A::Zeros)  = "Zero Operator"

zeros(B::AffineOperator, args...) = NestedLinearOperator(zeros, B, args...)


