import Base: +, -

type LinearOpSum{D3} <: LinearOp{D3}
	x::Array{OptVar}
	A::AbstractArray{LinearOp}
	mid::AbstractArray
	isTranspose::Bool
end

function size(A::LinearOpSum)
	dim1,dim2 = tuple(size.(A.x)...), size(A.A[1],2) 
	A.isTranspose ? (dim2,dim1) : (dim1,dim2)
end

fun_name(S::LinearOpSum) = "Sum of Linear Operators"

immutable SumSameVar{D1,D2} <: LinearOp{D1,D2}
	x::OptVar
	A::LinearOp
	B::LinearOp
	mid::AbstractArray
end
size(A::SumSameVar) = size(A.A)

fun_name(S::SumSameVar) = ((typeof(S.A) <: SumSameVar) == false ) ? 
fun_name(S.A)*" + "*fun_name(S.B) : "Sum of Linear Operators"

+(A::LinearOp, x::OptVar  ) = A+eye(x)
+(x::OptVar,   A::LinearOp) = eye(x)+A
+(x::OptVar,   y::OptVar  ) = eye(x)+eye(y)

-(A::LinearOp, x::OptVar  ) = A+(-eye(x))
-(x::OptVar,   A::LinearOp) = eye(x)+(-A)
-(x::OptVar,   y::OptVar  ) = eye(x)+(-eye(y))

function +{D1,D2,D3}(A::LinearOp{D1,D3}, B::LinearOp{D2,D3}) 
	if A.x == B.x
		if size(A) != size(B) DimensionMismatch() end
		mid = Array{D2}(size(B,2))
		return SumSameVar{D1,D2}(A.x,A,B,mid)
	else
		if size(A,2) != size(B,2) DimensionMismatch("operators must go to same space!") end
		mid = Array{D2}(size(B,2))
		LinearOpSum{D3}([A.x,B.x], [A,B], mid, false)
	end
end

function -{D1,D2,D3}(A::LinearOp{D1,D3}, B::LinearOp{D2,D3}) 
	return A+(-B)
end

function +{D1,D3}(A::LinearOpSum{D3},B::LinearOp{D1,D3})
	if A.isTranspose error("cannot sum to transpose") end
	if any(A.x .== B.x)
		for i = 1:length(A.A)
			if A.x[i] == B.x
				A.A[i] = A.A[i] + B 
			end
		end
	else
		push!(A.x,B.x)
		push!(A.A,B)
	end
	return A
end

function -{D1,D3}(A::LinearOpSum{D3},B::LinearOp{D1,D3})
	return A+(-B)
end

transpose{D1,D2}(S::SumSameVar{D1,D2}) = S.A'+ S.B' 
transpose{D2}(A::LinearOpSum{D2}) = LinearOpSum{D2}(A.x, A.A.', A.mid, A.isTranspose == false)

function A_mul_B!(y::AbstractArray,S::SumSameVar,b::AbstractArray) 
	A_mul_B!(S.mid,S.A,b)
	A_mul_B!(y  ,S.A,b)
	y .= (+).(S.mid,y)
end

function *{D2,T1<:AbstractArray}(A::LinearOpSum{D2},b::Array{T1,1}) 
	y = Array{D2}(size(A,2))
	A_mul_B!(y,A,b)
	return y
end

function *{D2}(A::LinearOpSum{D2},b::AbstractArray) 
	y = Array{AbstractArray,1}(length(A.A))
	for i = 1:length(A.A)
		y[i] = create_out(A.A[i]) 
	end
	A_mul_B!(y,A,b)
	return y
end

function A_mul_B!{T1<:AbstractArray}(y::AbstractArray,S::LinearOpSum,b::Array{T1,1}) 
	A_mul_B!(y,S.A[1],b[1])
	for i = 2:length(S.A)
		A_mul_B!(S.mid,S.A[i],b[i])
		y .= (+).(S.mid,y)
	end
end

create_out{D1,D2}(A::LinearOp{D1,D2}) = Array{D2}(size(A,2))

function A_mul_B!{T1<:AbstractArray}(y::Array{T1,1},S::LinearOpSum,b::AbstractArray) 
	for i = 1:length(S.A)
		A_mul_B!(y[i],S.A[i],b)
	end
end

function fun_dom{D3<:Real}(A::LinearOpSum{D3}) 
	str = ""
	for a in A.A str *= fun_D1(a) end
		str *= "→  ℝ^$(size(A,2))"
end

function fun_dom{D3<:Complex}(A::LinearOpSum{D3}) 
	str = ""
	for a in A.A str *= fun_D1(a) end
		str *= "→  ℂ^$(size(A,2))"
end

fun_D1{D1<:Real, D2}(A::LinearOp{D1,D2})    =  "ℝ^$(size(A,1)) "
fun_D1{D1<:Complex, D2}(A::LinearOp{D1,D2}) =  "ℂ^$(size(A,1)) "





