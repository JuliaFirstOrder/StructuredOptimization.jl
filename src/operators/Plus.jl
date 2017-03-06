import Base: +, -

immutable SumSameVar{D1,D2} <: LinearOp{D1,D2}
	x::OptVar
	A::LinearOp
	B::LinearOp
	mid::AbstractArray
	sign::Bool
end
size(A::SumSameVar) = size(A.A)

function SumSameVar{D1,D2}(x::OptVar,A::LinearOp{D1,D2},B::LinearOp{D1,D2}, sign::Bool) 
	mid = Array{D2}(size(A,2))
	return SumSameVar{D1,D2}(A.x,A,B,mid,sign)
end

fun_name(S::SumSameVar) = ((typeof(S.A) <: SumSameVar) == false ) ? 
fun_name(S.A)*(S.sign ? " + " : " - ")*fun_name(S.B) : "Sum of Linear Operators"

transpose{D1,D2}(S::SumSameVar{D1,D2}) = SumSameVar(S.x,S.A',S.B',S.sign)

function A_mul_B!(y::AbstractArray,S::SumSameVar,b::AbstractArray) 
	A_mul_B!(S.mid,S.A,b)
	A_mul_B!(y    ,S.B,b)
	S.sign ? y .= (+).(S.mid,y) : y .= (-).(S.mid,y)
end

+(A::LinearOp, x::OptVar  ) = A+eye(x)
+(x::OptVar,   A::LinearOp) = eye(x)+A
+(x::OptVar,   y::OptVar  ) = eye(x)+eye(y)
+(A::LinearOp) = A

-(A::LinearOp, x::OptVar  ) = A-eye(x)
-(x::OptVar,   A::LinearOp) = eye(x)-A
-(x::OptVar,   y::OptVar  ) = eye(x)-eye(y)
-{D1,D2}(A::LinearOp{D1,D2}) = Empty{D1,D2}(A.x)-A #trick to accept -A


function unsigned_sum{D1,D2,D3}(A::LinearOp{D1,D3}, B::LinearOp{D2,D3}, sign::Bool) 
	if A.x == B.x
		if size(A) == size(B) 
			return SumSameVar(A.x,A,B,sign)
		else
			return vcat(A,B,[true; sign])
		end
	else
		return hcat(A,B, [true; sign])
	end
end

+{D1,D2,D3}(A::LinearOp{D1,D3}, B::LinearOp{D2,D3}) = unsigned_sum(A,B,true ) 
-{D1,D2,D3}(A::LinearOp{D1,D3}, B::LinearOp{D2,D3}) = unsigned_sum(A,B,false)

function unsigned_sum{D1,D3}(A::HCAT{D3},B::LinearOp{D1,D3}, sign::Bool)
	if size(A,2) != size(B,2) DimensionMismatch("Operators must share codomain") end
	if any(A.x .== B.x)
		for i = 1:length(A.A)
			if A.x[i] == B.x
				sign ? A.A[i] = A.A[i] + B : A.A[i] = A.A[i] - B
				break
			end
		end
	else
		push!(A.x,B.x)
		push!(A.A,B)
		push!(A.sign,sign)
	end
	return A
end

+{D1,D3}(A::HCAT{D3},B::LinearOp{D1,D3}) = unsigned_sum(A,B,true )
-{D1,D3}(A::HCAT{D3},B::LinearOp{D1,D3}) = unsigned_sum(A,B,false)








