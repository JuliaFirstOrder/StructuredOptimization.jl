import Base: +, -
export operator, adjoint, variable

type Affine <: AffineOperator
	x::Vector{OptVar}
	A::LinearOperator
	At::LinearOperator
	b::Nullable{Vector{AbstractArray}}
end
variable(A::Affine) = A.x
operator(A::Affine) = A.A
adjoint(A::Affine)  = A.At

  domainType(A::AffineOperator) =   domainType(operator(A))
codomainType(A::AffineOperator) = codomainType(operator(A))

fun_name(A::Affine) = "Affine "*fun_name(A.A)
fun_dom(A::Affine)  = fun_dom(A.A)

-(A::LinearOperator) = emptyop(A)-A #trick to accept -A
-(A::Affine) = (A.A = -A.A; A.At = -A.At; return A) 
-(x::OptVar) = -(eye(x)) 

+(A::Affine,b::AbstractArray)  =  
(isnull(A.b) ? A.b = Nullable([b]) : A.b = Nullable(get(A.b)+[b]); return A)
+(b::AbstractArray,A::Affine)  = A+b 
-(A::Affine,b::AbstractArray)  = A+(-b)
-(b::AbstractArray, A::Affine) =  (-A)+b

+(x::OptVar,b::AbstractArray) = eye(x)+b 
-(x::OptVar,b::AbstractArray) = eye(x)-b
+(b::AbstractArray,x::OptVar) = b+eye(x) 
-(b::AbstractArray,x::OptVar) = b+(-eye(x))

function evaluate!(y::AbstractArray, A::Affine, x::AbstractArray) 
	A_mul_B!(y,A.A,x)
	isnull(A.b) ? nothing : y .+= get(A.b)[1]
end

function (A::Affine)(x::AbstractArray) 
	y = Array{codomainType(operator(A))}(size(operator(A),2))
	evaluate!(y,A,x)
	return y
end

function Base.show{Op <: AffineOperator }(io::IO, f::Op)
  println(io, "description : ", fun_name(f))
  println(io, "domain      : ", fun_dom(f))
end
