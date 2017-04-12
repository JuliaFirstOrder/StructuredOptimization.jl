import Base: +,-

immutable TiltedAffine <: AffineOperator
	A::Affine
	b::AbstractArray

	function TiltedAffine(A,b)
		size(b)!= size(A,1) && throw(DimensionMismatch()) 
		codomainType(A) != eltype(b) && throw(DomainError())
		new(A,b)
	end
end

variable(A::TiltedAffine) = variable(A.A)
operator(A::TiltedAffine) = operator(A.A)
adjoint( A::TiltedAffine) = adjoint( A.A)
tilt(    A::TiltedAffine) = A.b 

domainType(A::TiltedAffine) =   domainType(A::Affine)
codomainType(A::TiltedAffine) = codomainType(A::Affine)
size(A::TiltedAffine, args...) = size(A.A, args...)

# Constructors

(+)(A::Affine,b::AbstractArray)       = TiltedAffine( A, b)
(-)(A::Affine,b::AbstractArray)       = TiltedAffine( A,-b)
(+)(b::AbstractArray,A::Affine)       = TiltedAffine( A, b)
(-)(b::AbstractArray,A::Affine)       = TiltedAffine(-A, b)
(+)(A::TiltedAffine,b::AbstractArray) = TiltedAffine( A.A, A.b+b)
(-)(A::TiltedAffine,b::AbstractArray) = TiltedAffine( A.A, A.b-b)
(+)(b::AbstractArray,A::TiltedAffine) = TiltedAffine( A.A, A.b+b)
(-)(b::AbstractArray,A::TiltedAffine) = TiltedAffine(-A.A, b-A.b)

(+)(A::TiltedAffine,B::Affine) = TiltedAffine(A.A+B, A.b )
(-)(A::TiltedAffine,B::Affine) = TiltedAffine(A.A-B, A.b )
(+)(A::Affine,B::TiltedAffine) = TiltedAffine(A+B.A, B.b )
(-)(A::Affine,B::TiltedAffine) = TiltedAffine(A-B.A,-B.b )

(+)(A::TiltedAffine,B::TiltedAffine) = TiltedAffine(A.A+B.A, A.b+B.b )
(-)(A::TiltedAffine,B::TiltedAffine) = TiltedAffine(A.A-B.A, A.b-B.b )

(+)(x::Variable,b::AbstractArray) = eye(x)+b
(-)(x::Variable,b::AbstractArray) = eye(x)-b
(+)(b::AbstractArray,x::Variable) = b+eye(x)
(-)(b::AbstractArray,x::Variable) = b-eye(x)

# Operators

function evaluate!(y::AbstractArray, A::TiltedAffine, x::AbstractArray)
	evaluate!(y, A.A, x)
	y .+= A.b
end

function (A::TiltedAffine)(x::AbstractArray)
	y = A.A(x)
	y .+= A.b
	return y
end

#sorting operators

sort_and_expand(x, A::TiltedAffine) = TiltedAffine(sort_and_expand(x, A.A),A.b) 

# Printing stuff

fun_name(A::TiltedAffine) = fun_name(A.A)*" + b"
fun_type(A::TiltedAffine) = fun_type(A.A)

function Base.show(io::IO, f::TiltedAffine)
  println(io, "description : Tilted Affine")
  println(io, "operator    : ", fun_name(f) )
  println(io, "type        : ", fun_type(f))
end
