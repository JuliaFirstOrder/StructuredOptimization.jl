import Base: +, -

immutable SumSameVar{D1,D2} <: LinearOperator{D1,D2}
	A::Vector{LinearOperator}
	mid::AbstractArray
	sign::Vector{Bool}
end
size(S::SumSameVar) = size(S.A[1])

function unsigned_sum{D1,D2}(A::LinearOperator{D1,D2},B::LinearOperator{D1,D2}, sign::Bool) 
	if size(A) != size(B) dimension_error(A,B) end
	mid = Array{D2}(size(A,2))
	return SumSameVar{D1,D2}([A,B],mid,[true,sign])
end

function unsigned_sum{D1,D2}(S::SumSameVar{D1,D2},B::LinearOperator{D1,D2}, sign::Bool) 
	if size(S) != size(B) dimension_error(S,B) end
	A = copy(S.A)
	s = copy(S.sign)
	push!(A,B)
	push!(s,sign)
	return SumSameVar{D1,D2}(A,S.mid,s)
end

+{D1,D2}(A::LinearOperator{D1,D2}, B::LinearOperator{D1,D2}) = unsigned_sum(A,B,true ) 
-{D1,D2}(A::LinearOperator{D1,D2}, B::LinearOperator{D1,D2}) = unsigned_sum(A,B,false)

function fun_name(S::SumSameVar) 
	if length(S.A) == 2  
		fun_name(S.A[1])*(S.sign[2] ? " + " : " - ")*fun_name(S.A[2]) 
	else 
		"Sum of Linear Operators"
	end
end

transpose{D1,D2}(S::SumSameVar{D1,D2}) = SumSameVar{D2,D1}((S.A.')[:],Array{D1}(size(S,1)),S.sign)

function A_mul_B!(y::AbstractArray,S::SumSameVar,b::AbstractArray) 
	A_mul_B!(y,S.A[1],b)
	for i = 2:length(S.A)
		A_mul_B!(S.mid,S.A[i],b)
		S.sign[i] ? y .= (+).(y,S.mid) : y .= (-).(y,S.mid)
	end
end

-(A::LinearOperator) = zeros(A)-A #trick to accept -A
-(A::Affine) = isnull(A.b) ? Affine(A.x,-A.A,-A.At,A.b) : Affine(A.x,-A.A,-A.At,Nullable(-get(A.b)))
-(x::OptVar) = -(eye(x)) 

+(A::Affine,b::AbstractArray)  =  
(isnull(A.b) ? Affine(A.x,A.A,A.At,Nullable([b])) : Affine(A.x,A.A,A.At,Nullable(get(A.b)+[b])))
+(b::AbstractArray,A::Affine)  = A+b 
-(A::Affine,b::AbstractArray)  = A+(-b)
-(b::AbstractArray, A::Affine) =  (-A)+b

+(x::OptVar,b::AbstractArray) = eye(x)+b 
-(x::OptVar,b::AbstractArray) = eye(x)-b
+(b::AbstractArray,x::OptVar) = b+eye(x) 
-(b::AbstractArray,x::OptVar) = b+(-eye(x))

+(x::OptVar,A::Affine) = eye(x)+A 
-(x::OptVar,A::Affine) = eye(x)-A
+(A::Affine,x::OptVar) = A+eye(x) 
-(A::Affine,x::OptVar) = A-eye(x)

function unsigned_sum(A::Affine,B::Affine, sign::Bool)  
	if isnull(A.b) == true && isnull(B.b) == true
		b = Nullable{Vector{AbstractArray}}()
	elseif isnull(A.b) == false && isnull(B.b) == true
		b = Nullable(get(A.b))
	elseif isnull(A.b) == true && isnull(B.b) == false
		sign ? b = Nullable(get(B.b)) : b = Nullable(-get(B.b))
	elseif isnull(A.b) == false && isnull(B.b) == false
		sign ? b = Nullable(get(A.b)+get(B.b)) : b = Nullable(get(A.b)-get(B.b))
	end

	if variable(A) == variable(B)
		if size(operator(A)) == size(operator(B)) 
			S = unsigned_sum(operator(A),operator(B),sign)
			x = A.x
		else
			dimension_error(operator(A),operator(B))
		end
	else
		if size(operator(A),2) == size(operator(B),2) 
			S,x = unsigned_sum(variable(A),operator(A),variable(B),operator(B),sign) 
		else
			dimension_error(operator(A),operator(B))
		end
	end

	return Affine(x, S, S',b)
end

+(A::Affine, B::Affine) = unsigned_sum(A,B,true ) 
-(A::Affine, B::Affine) = unsigned_sum(A,B,false)

function unsigned_sum{D1,D2,D3}(xa::Vector{AbstractOptVar}, A::LinearOperator{D1,D3}, 
				xb::Vector{AbstractOptVar}, B::LinearOperator{D2,D3}, sign::Bool)
	return hcat(A,B,[true,sign]), [xa[1],xb[1]]

end

function unsigned_sum{D1,D2,
		      T1<:AbstractOptVar,
		      T2<:AbstractOptVar}(xa::Vector{T1}, A::HCAT{D2}, 
			    		  xb::Vector{T2}, B::LinearOperator{D1,D2}, sign::Bool)
	H = copy(A.A)
	x = copy(xa)
	s = copy(A.sign)

	if any(x .== xb[1])
		idx = find(x .== xb[1])[1]
		if s[idx] == true
			sign ? H[idx] = H[idx] + B : H[idx] =  H[idx] - B
		else
			sign ? H[idx] = H[idx] - B : H[idx] =  H[idx] + B
		end
	else
		push!(x,xb[1])
		push!(H,B)
		push!(s,sign)
	end
	return HCAT{D2}(H,A.mid,s), x
end

unsigned_sum{D1,D2}(xa::Vector{AbstractOptVar}, A::LinearOperator{D1,D2}, 
		    xb::Vector{AbstractOptVar}, B::HCAT{D2}, sign::Bool) = unsigned_sum(xb,B,xa,A,sign)

function unsigned_sum{D2,
		      T1<:AbstractOptVar,
		      T2<:AbstractOptVar}(xa::Vector{T1}, A::HCAT{D2}, 
			    		  xb::Vector{T2}, B::HCAT{D2}, sign::Bool)

	H, x = unsigned_sum(xa,A,[xb[1]],B.A[1],sign)
	for i = 2:length(xb)
		H, x = unsigned_sum(x,H,[xb[i]],B.A[i],sign)
	end
	return H, x
end

dimension_error(A::LinearOperator,B::LinearOperator) =		
throw(DimensionMismatch("cannot sum operator of size $(size(A)) with operator of size$(size(B))"))





