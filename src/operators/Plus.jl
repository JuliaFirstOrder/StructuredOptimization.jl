import Base: +, -

immutable SumSameVar{D1,D2} <: LinearOp{D1,D2}
	x::OptVar
	A::LinearOp
	B::LinearOp
	mid::AbstractArray
	sign::Bool
end
size(A::SumSameVar) = size(A.A)

fun_name(S::SumSameVar) = ((typeof(S.A) <: SumSameVar) == false ) ? 
fun_name(S.A)*(S.sign ? " + " : " - ")*fun_name(S.B) : "Sum of Linear Operators"

+(A::LinearOp, x::OptVar  ) = A+eye(x)
+(x::OptVar,   A::LinearOp) = eye(x)+A
+(x::OptVar,   y::OptVar  ) = eye(x)+eye(y)
+(A::LinearOp) = A

-(A::LinearOp, x::OptVar  ) = A-eye(x)
-(x::OptVar,   A::LinearOp) = eye(x)-A
-(x::OptVar,   y::OptVar  ) = eye(x)-eye(y)
-{D1,D2}(A::LinearOp{D1,D2}) = Empty{D1,D2}(A.x)-A #trick to accept -A

function +{D1,D2,D3}(A::LinearOp{D1,D3}, B::LinearOp{D2,D3}) 
	if A.x == B.x
		if size(A) != size(B) DimensionMismatch() end
		mid = Array{D3}(size(B,2))
		return SumSameVar{D1,D3}(A.x,A,B,mid,true)
	else
		hcat(A,B, [true; true])
	end
end

function -{D1,D2,D3}(A::LinearOp{D1,D3}, B::LinearOp{D2,D3}) 
	if A.x == B.x
		if size(A) != size(B) DimensionMismatch() end
		mid = Array{D3}(size(B,2))
		return SumSameVar{D1,D3}(A.x,A,B,mid,false)
	else
		hcat(A,B, [true; false])
	end
end

function +{D1,D3}(A::HCAT{D3},B::LinearOp{D1,D3})
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
		push!(A.sign,true)
	end
	return A
end

function -{D1,D3}(A::HCAT{D3},B::LinearOp{D1,D3})
	if A.isTranspose error("cannot sum to transpose") end
	if any(A.x .== B.x)
		for i = 1:length(A.A)
			if A.x[i] == B.x
				A.A[i] = A.A[i] - B 
			end
		end
	else
		push!(A.x,B.x)
		push!(A.A,B)
		push!(A.sign,false)
	end
	return A
end

transpose{D1,D2}(S::SumSameVar{D1,D2}) = S.sign ? S.A'+ S.B' : S.A'- S.B'

function A_mul_B!(y::AbstractArray,S::SumSameVar,b::AbstractArray) 
	A_mul_B!(S.mid,S.A,b)
	A_mul_B!(y    ,S.B,b)
	S.sign ? y .= (+).(S.mid,y) : y .= (-).(S.mid,y)
end








