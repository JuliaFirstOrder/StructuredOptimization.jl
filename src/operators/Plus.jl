import Base: +, -

immutable SumSameVar{D1,D2} <: LinearOperator{D1,D2}
	A::LinearOperator
	B::LinearOperator
	mid::AbstractArray
	sign::Bool
end
size(A::SumSameVar) = size(A.A)

function SumSameVar{D1,D2}(A::LinearOperator{D1,D2},B::LinearOperator{D1,D2}, sign::Bool) 
	mid = Array{D2}(size(A,2))
	return SumSameVar{D1,D2}(A,B,mid,sign)
end

fun_name(S::SumSameVar) = ((typeof(S.A) <: SumSameVar) == false ) ? 
fun_name(S.A)*(S.sign ? " + " : " - ")*fun_name(S.B) : "Sum of Linear Operators"

transpose{D1,D2}(S::SumSameVar{D1,D2}) = SumSameVar(S.A',S.B',S.sign)

function A_mul_B!(y::AbstractArray,S::SumSameVar,b::AbstractArray) 
	A_mul_B!(S.mid,S.A,b)
	A_mul_B!(y    ,S.B,b)
	S.sign ? y .= (+).(S.mid,y) : y .= (-).(S.mid,y)
end


function unsigned_sum{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3}, sign::Bool) 
		
	if size(A) == size(B) 
		return SumSameVar(A,B,sign)
	else
		DimensionMismatch("cannot sum operator of size $(size(A)) with operator of size$(size(B))")
	
	end
end

+{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3}) = unsigned_sum(A,B,true ) 
-{D1,D2,D3}(A::LinearOperator{D1,D3}, B::LinearOperator{D2,D3}) = unsigned_sum(A,B,false)

function unsigned_sum{D1,D3}(A::HCAT{D3},B::LinearOperator{D1,D3}, sign::Bool)
	if size(A,2) != size(B,2) DimensionMismatch("Operators must share codomain") end
	if any(variable(A) .== variable(B))
		for i = 1:length(A.A)
			if variable(A)[i] == variable(B)
				sign ? A.A[i] = A.A[i] + B : A.A[i] = A.A[i] - B
				break
			end
		end
	else
		push!(A.A,B)
		push!(A.sign,sign)
	end
	return A
end

+{D1,D3}(A::HCAT{D3},B::LinearOperator{D1,D3}) = unsigned_sum(A,B,true )
-{D1,D3}(A::HCAT{D3},B::LinearOperator{D1,D3}) = unsigned_sum(A,B,false)








