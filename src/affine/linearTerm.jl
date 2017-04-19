import Base: +, -
export operator, adjoint, variable, LinearTerm

immutable LinearTerm <: AbstractAffineTerm
	x::Vector{AbstractVariable}
	L::LinearOperator
	Lt::LinearOperator
	function LinearTerm(x,L)
		(length(x) == 1 ? size(x[1]) != size(L,2) :
                 any(size.(x) != [size(L,2)...]) ) && throw(DimensionMismatch())
		new(x,L,L')
	end
end

variable(A::LinearTerm) = A.x
operator(A::LinearTerm) = A.L
adjoint(A::LinearTerm)  = A.Lt
tilt(A::LinearTerm)     =  0.

domainType(A::LinearTerm) =   domainType(operator(A))
codomainType(A::LinearTerm) = codomainType(operator(A))
size(A::LinearTerm, args...) = size(operator(A), args...)

# Constructors

LinearTerm(x::Variable,L::LinearOperator) = LinearTerm([x],L)

+(A::LinearTerm) = A
-(A::LinearTerm) = LinearTerm(variable(A),-operator(A))

+(A::LinearTerm, B::LinearTerm) =  all(variable(A) == variable(B)) ? LinearTerm(variable(A), A.L+B.L) :
LinearTerm(unsigned_sum(variable(A),operator(A),variable(B),operator(B), true)...)
#see utils down here for unsigned_sum
# other constructors in AffineConstructors and ComposeAffine

-(A::LinearTerm, B::LinearTerm) =  all(variable(A) == variable(B)) ? LinearTerm(variable(A), A.L-B.L) :
LinearTerm(unsigned_sum(variable(A),operator(A),variable(B),operator(B), false)...)

-(x::Variable) = -eye(x)
+(x::Variable,A::LinearTerm) = eye(x)+A
-(x::Variable,A::LinearTerm) = eye(x)-A
+(A::LinearTerm,x::Variable) = A+eye(x)
-(A::LinearTerm,x::Variable) = A-eye(x)

+(x::Variable,y::Variable) = eye(x)+eye(y)
-(x::Variable,y::Variable) = eye(x)-eye(y)

# Operators

function evaluate!(y::AbstractArray, A::LinearTerm, x::AbstractArray)
	A_mul_B!(y,operator(A),x)
end

function (A::LinearTerm)(x::AbstractArray)
	y = operator(A)*x
	return y
end

#sorting operators

function sort_and_expand{T<:AbstractVariable}(x::Vector{T}, A::LinearTerm)
	if all(x == A.x)
		return A
	else
		dim = size(operator(A),1)

		H = Vector{LinearOperator}(length(x))
		[H[i] = Zeros(dim,size(x[i]))  for i in eachindex(H) ]
		for i in eachindex(x)
			if any(A.x .== x[i])
				idx = find(A.x .== x[i])[1]
				H[i] = extract_operator(operator(A),idx)
			end
		end
		H = hcat(H...)
		return LinearTerm(x,H)
	end
end


extract_operator(L::LinearOperator, idx::Int64) = L
extract_operator(L::HCAT, idx::Int64) = L.A[idx]

# Printing stuff

fun_name(A::LinearTerm) = fun_name(operator(A))
fun_type(A::LinearTerm) = fun_type(operator(A))

function Base.show(io::IO, f::LinearTerm)
  println(io, "description : LinearTerm")
  println(io, "operator    : ", fun_name(operator(f)))
  print(  io, "type        : ", fun_type(operator(f)))
end

# utils

function unsigned_sum(xa::Vector{AbstractVariable}, A::LinearOperator,
		      xb::Vector{AbstractVariable}, B::LinearOperator, sign::Bool)
	sign ? ([xa[1],xb[1]], hcat(A,B)) : ([xa[1],xb[1]], hcat(A,-B))
end

function unsigned_sum(xa::Vector{AbstractVariable}, A::HCAT,
		      xb::Vector{AbstractVariable}, B::LinearOperator, sign::Bool)
	H = copy(A.A)
	x = copy(xa)

	if any(x .== xb[1])
		idx = find(x .== xb[1])[1]
		sign ? H[idx] = H[idx] + B : H[idx] =  H[idx] - B
	else
		push!(x,xb[1])
		sign ? push!(H,B) : push!(H,-B)
	end
	return x, HCAT(H,A.mid)
end

unsigned_sum(xa::Vector{AbstractVariable}, A::LinearOperator,
	     xb::Vector{AbstractVariable}, B::HCAT, sign::Bool) = unsigned_sum(xb,B,xa,A,sign)

function unsigned_sum(xa::Vector{AbstractVariable}, A::HCAT,
		      xb::Vector{AbstractVariable}, B::HCAT, sign::Bool)

	x, H = unsigned_sum(xa,A,[xb[1]],B.A[1],sign)
	for i = 2:length(xb)
		x, H = unsigned_sum(x,H,[xb[i]],B.A[i],sign)
	end
	return x,H
end
