immutable NestedLinearOperator{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	A::Vector{LinearOperator}
	mid::Vector{AbstractArray}       # memory in the middle of the operators
	NestedLinearOperator(sign,A,mid) = new(sign,A,mid)
	NestedLinearOperator(     A,mid) = new(true,A,mid)
end

size(A::NestedLinearOperator) = (size(A.A[1],1), size(A.A[end],2))
-{D1,D2}(A::NestedLinearOperator{D1,D2}) =
NestedLinearOperator{D1,D2}(false == sign(A), A.A, A.mid)

function NestedLinearOperator(f::Function,B::AffineOperator, args...)
	mid = Array{codomainType(operator(B))}(size(operator(B),2))
	(f == *) ? A = f(args[1], Variable(mid)) : A = f(Variable(mid), args...)
	N = NestedLinearOperator(operator(A),operator(B),mid)
	b = Nullable{AbstractArray}()
	isnull(B.b) ? nothing : b = adjoint(A)*get(B.b)
	Affine(variable(B),N,N',b)
end

*(A::LinearOperator, B::LinearOperator) = NestedLinearOperator(A,B)
*{T<:Number}(a::T, B::LinearOperator)   = NestedLinearOperator(DiagOp(a,size(B,2)...),B)
*{E<:IdentityOperator}(A::E, B::LinearOperator) = sign(A) ? B : -B

.*(A::LinearOperator, B::LinearOperator) = NestedLinearOperator(A,B)
.*{E<:IdentityOperator}(A::E, B::LinearOperator) = A*B

NestedLinearOperator{D1,Dm,D2}(A::LinearOperator{Dm,D2},B::LinearOperator{D1,Dm}) =
NestedLinearOperator(A,B,Array{Dm}(size(B,2)))


NestedLinearOperator{D1,Dm,D2}(A::LinearOperator{Dm,D2},
				B::LinearOperator{D1,Dm},mid::AbstractArray)  =
NestedLinearOperator{D1,D2}(sign(A) == sign(B), [B,A], [mid])

NestedLinearOperator{D1,Dm,D2}(A::NestedLinearOperator{Dm,D2},
				B::LinearOperator{D1,Dm},mid::AbstractArray)  =
NestedLinearOperator{D1,D2}(sign(A) == sign(B), [B,A.A...], [mid,A.mid...])

NestedLinearOperator{D1,Dm,D2}(A::LinearOperator{Dm,D2},
				B::NestedLinearOperator{D1,Dm},mid::AbstractArray)  =
NestedLinearOperator{D1,D2}(sign(A) == sign(B), [B.A...,A], [B.mid...,mid])

NestedLinearOperator{D1,Dm,D2}(A::NestedLinearOperator{Dm,D2},
				B::NestedLinearOperator{D1,Dm},mid::AbstractArray) =
NestedLinearOperator{D1,D2}(sign(A) == sign(B), [B.A...,A.A...], [B.mid...,mid,A.mid...])


function transpose{D1,D2}(N::NestedLinearOperator{D1,D2})
	NestedLinearOperator{D2,D1}(sign(N),flipdim((N.A.')[:],1),flipdim(N.mid,1))
end

function uA_mul_B!{D1,D2}(y::AbstractArray,N::NestedLinearOperator{D1,D2},b::AbstractArray)
	uA_mul_B!(N.mid[1],N.A[1],b)
	for i = 2:length(N.A)-1
		uA_mul_B!(N.mid[i],N.A[i], N.mid[i-1])
	end
	uA_mul_B!(y,N.A[length(N.A)], N.mid[length(N.A)-1])
end

fun_name(N::NestedLinearOperator) = length(N.A) == 2 ? fun_name(N.A[2])*" * "*fun_name(N.A[1]) : "Nested Linear Operator"
