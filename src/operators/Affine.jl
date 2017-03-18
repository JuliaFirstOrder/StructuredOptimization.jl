import Base: +, -

immutable Affine{D2} <: AffineOperator
	x::Vector{OptVar}
	A::LinearOperator
	b::Vector{Nullable{AbstractArray}}
end
variable(A::Affine) = variable(A.A)

fun_name(A::Affine) = "Affine "*fun_name(A.A)
fun_dom(A::Affine)  = fun_dom(A.A)

+(A::Affine,b::AbstractArray)  =  
(isnull(A.b[1]) ? A.b .= [Nullable( b)] : A.b .= [Nullable(get(A.b)+b)]; return A)
-(A::Affine,b::AbstractArray)  = 
(isnull(A.b[1]) ? A.b .= [Nullable(-b)] : A.b .= [Nullable(get(A.b)-b)]; return A)

+(x::OptVar,b::AbstractArray) = Affine([x],eye(x),[ b])
-(x::OptVar,b::AbstractArray) = Affine([x],eye(x),[-b])
+(b::AbstractArray,x::OptVar) = Affine([x],eye(x),[ b])
-(b::AbstractArray,x::OptVar) = Affine([x],-eye(x),[ b])

function evaluate!(y::AbstractArray, A::Affine, x::AbstractArray) 
	A_mul_B!(y,A.A,x)
	isnull(A.b[1]) ? nothing : y .+= get(A.b[1])
end

function (A::Affine{D2}){D2}(x::AbstractArray) 
	y = Array{D2}(size(A.A,2))
	evaluate!(y,A,x)
	return y
end

function Base.show{Op <: AffineOperator }(io::IO, f::Op)
  println(io, "description : ", fun_name(f))
  println(io, "domain      : ", fun_dom(f))
end
