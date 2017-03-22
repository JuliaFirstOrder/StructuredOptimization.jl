import Base: +, -
export operator, adjoint, variable, Affine

immutable Affine <: AffineOperator
	x::Vector{AbstractOptVar}
	A::LinearOperator
	At::LinearOperator
	b::Nullable{Vector{AbstractArray}}
end
variable(A::Affine) = A.x
operator(A::Affine) = A.A
adjoint(A::Affine)  = A.At

#TODO add checks of dimension of OptVar  
Affine(x::OptVar,A::LinearOperator) = Affine([x],A)
Affine{T<:AbstractOptVar}(x::Vector{T}, A::LinearOperator) = Affine(x, A, A',Nullable{Vector{AbstractArray}}()) 

  domainType(A::AffineOperator) =   domainType(operator(A))
codomainType(A::AffineOperator) = codomainType(operator(A))

fun_name(A::Affine) = "Affine "*fun_name(A.A)
fun_dom(A::Affine)  = fun_dom(A.A)

function evaluate!(y::AbstractArray, A::Affine, x::AbstractArray) 
	A_mul_B!(y,A.A,x)
	add_b!(y,A)
end

function (A::Affine)(x::AbstractArray) 
	y = operator(A)*x
	add_b!(y,A)
	return y
end

add_b!(y::AbstractArray,A::Affine) = isnull(A.b) ? nothing : y .+= get(A.b)[1]

function Base.show{Op <: AffineOperator }(io::IO, f::Op)
  println(io, "description : ", fun_name(f))
  println(io, "domain      : ", fun_dom(f))
end

#sorting operators

import Base: sort, issorted, sortperm

issorted(A::Affine) = issorted(A.x, by = object_id)
sortperm(A::Affine) = sortperm(A.x, by = object_id)

function sort(A::Affine)
	if length(A.x) == 1
		return A
	else
		p = sortperm(A.x, by = object_id)
		x  = A.x[p]
		H  = sort(A.A, p)
		return Affine(x,H,H',A.b)
	end
end

function sort_and_expand{T<:AbstractOptVar}(x::Vector{T}, A::Affine)
	if all(x == A.x)
		return A
	else
		D3   = codomainType(operator(A))
		dim2 = size(operator(A),2)

		H = Vector{LinearOperator}(length(x))
		[H[i] = Zeros(size(x[i]),dim2)  for i in eachindex(H) ]
		s = ones(Bool,length(x))
		for i in eachindex(x)
			if any(A.x .== x[i])
				idx = find(A.x .== x[i])[1]
				H[i],s[i] = extract_operator(operator(A),idx)
			end
		end
		H = hcat(H...)
		return Affine(x,H,H',A.b)
	end
end

extract_operator(A::LinearOperator, idx::Int64) = A,true
















