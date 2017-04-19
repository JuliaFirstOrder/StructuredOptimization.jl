import Base: +,-

immutable AffineTerm <: AbstractAffineTerm
	A::LinearTerm
	b::AbstractArray

	function AffineTerm(A,b)
		size(b)!= size(A,1) && throw(DimensionMismatch())
		codomainType(A) != eltype(b) && throw(DomainError())
		new(A,b)
	end
end

variable(A::AffineTerm) = variable(A.A)
operator(A::AffineTerm) = operator(A.A)
adjoint( A::AffineTerm) = adjoint( A.A)
tilt(    A::AffineTerm) = A.b

domainType(A::AffineTerm) =   domainType(A::LinearTerm)
codomainType(A::AffineTerm) = codomainType(A::LinearTerm)
size(A::AffineTerm, args...) = size(A.A, args...)

# Constructors

(+)(A::LinearTerm,b::AbstractArray)       = AffineTerm( A, b)
(-)(A::LinearTerm,b::AbstractArray)       = AffineTerm( A,-b)
(+)(b::AbstractArray,A::LinearTerm)       = AffineTerm( A, b)
(-)(b::AbstractArray,A::LinearTerm)       = AffineTerm(-A, b)
(+)(A::AffineTerm,b::AbstractArray) = AffineTerm( A.A, A.b+b)
(-)(A::AffineTerm,b::AbstractArray) = AffineTerm( A.A, A.b-b)
(+)(b::AbstractArray,A::AffineTerm) = AffineTerm( A.A, A.b+b)
(-)(b::AbstractArray,A::AffineTerm) = AffineTerm(-A.A, b-A.b)

(+)(A::AffineTerm,B::LinearTerm) = AffineTerm(A.A+B, A.b )
(-)(A::AffineTerm,B::LinearTerm) = AffineTerm(A.A-B, A.b )
(+)(A::LinearTerm,B::AffineTerm) = AffineTerm(A+B.A, B.b )
(-)(A::LinearTerm,B::AffineTerm) = AffineTerm(A-B.A,-B.b )

(+)(A::AffineTerm,B::AffineTerm) = AffineTerm(A.A+B.A, A.b+B.b )
(-)(A::AffineTerm,B::AffineTerm) = AffineTerm(A.A-B.A, A.b-B.b )

(+)(x::Variable,b::AbstractArray) = eye(x)+b
(-)(x::Variable,b::AbstractArray) = eye(x)-b
(+)(b::AbstractArray,x::Variable) = b+eye(x)
(-)(b::AbstractArray,x::Variable) = b-eye(x)

#special cases

(+)(x::Variable,b::Number) = b == 0 ? eye(x) : error("cannot sum $(typeof(x)) with $(typeof(b))")
(-)(x::Variable,b::Number) = b == 0 ? eye(x) : error("cannot sum $(typeof(x)) with $(typeof(b))")
(+)(A::AbstractAffineTerm,b::Number) = b == 0 ? A : error("cannot sum $(typeof(A)) with $(typeof(b))")
(-)(A::AbstractAffineTerm,b::Number) = b == 0 ? A : error("cannot sum $(typeof(A)) with $(typeof(b))")

# Operators

function evaluate!(y::AbstractArray, A::AffineTerm, x::AbstractArray)
	evaluate!(y, A.A, x)
	y .+= A.b
end

function (A::AffineTerm)(x::AbstractArray)
	y = A.A(x)
	y .+= A.b
	return y
end

#sorting operators

sort_and_expand(x, A::AffineTerm) = AffineTerm(sort_and_expand(x, A.A),A.b)

# Printing stuff

fun_name(A::AffineTerm) = fun_name(A.A)*" + b"
fun_type(A::AffineTerm) = fun_type(A.A)

function Base.show(io::IO, f::AffineTerm)
  println(io, "description : Tilted LinearTerm")
  println(io, "operator    : ", fun_name(f) )
  print(  io, "type        : ", fun_type(f))
end
