import Base: +, -

"""
`+(A::AbstractOperator, ex::AbstractExpression)`

Add two expressions. 

# Examples

```julia
julia> x,y = Variable(5), Variable(5)
(Variable(Float64, (5,)), Variable(Float64, (5,)))

julia> ex1 = x+y

julia> z = Variable(2)
Variable(Float64, (2,))

julia> ex2 = randn(5,2)*z

```

Notice that in order for two expressions to be added toghether their associate `AbstractOperator` 
must have the same codomain:

```julia
julia> operator(ex1)
[I,I]  ℝ^5  ℝ^5 -> ℝ^5 

julia> operator(ex2)
▒  ℝ^2 -> ℝ^5 

julia> ex3 = ex1 + ex2

```

It is also possible to use broadcasted addition:
```julia 
julia> z = Variable(1)
Variable(Float64, (1,))

julia> ex3.+z

```

"""
function (+)(a::AbstractExpression, b::AbstractExpression)
	A = convert(Expression,a)
	B = convert(Expression,b)
	d = displacement(A)+displacement(B)
	if variables(A) == variables(B)
		return Expression{length(A.x)}(A.x,operator(A)+operator(B),d)
	else
		opA = operator(A)
		xA = variables(A)
		opB = operator(B)
		xB = variables(B)

		xNew, opNew = Usum_op(xA,xB,opA,opB,true)
		return Expression{length(xNew)}(xNew,opNew,d)
	end

end
# sum expressions

function (-)(a::AbstractExpression, b::AbstractExpression)
	A = convert(Expression,a)
	B = convert(Expression,b)
	d = displacement(A)-displacement(B)
	if variables(A) == variables(B)
		return Expression{length(A.x)}(A.x,operator(A)-operator(B),d)
	else
		opA = operator(A)
		xA = variables(A)
		opB = operator(B)
		xB = variables(B)

		xNew, opNew = Usum_op(xA,xB,opA,opB,false)
		return Expression{length(xNew)}(xNew,opNew,d)
	end

end

#unsigned sum operators with single variables
function Usum_op(xA::Tuple{Variable},
		 xB::Tuple{Variable},
		 A::AbstractOperator,
		 B::AbstractOperator,sign::Bool) 
	xNew  = (xA...,xB...)
	opNew = sign ? hcat(A,B) : hcat(A,-B)
	return xNew, opNew
end

#unsigned sum: HCAT + AbstractOperator
function Usum_op(xA::NTuple{N,Variable},
		 xB::Tuple{Variable},
		 A::L1,
		 B::AbstractOperator,sign::Bool) where {N, M, L1<:HCAT{M,N}}
	if xB[1] in xA
		idx = findfirst(xA.==xB[1])
		S = sign ? A[idx]+B : A[idx]-B
		xNew = xA
		opNew = hcat(A[1:idx-1],S,A[idx+1:N]  )
	else
		xNew  = (xA...,xB...)
		opNew = sign ? hcat(A,B) : hcat(A,-B)
	end

	return xNew, opNew
end

#unsigned sum: AbstractOperator+HCAT
function Usum_op(xA::Tuple{Variable},
		 xB::NTuple{N,Variable},
		 A::AbstractOperator,
		 B::L2,sign::Bool) where {N, M, L2<:HCAT{M,N}}
	if xA[1] in xB
		idx = findfirst(xA.==xB[1])
		S = sign ? A+B[idx] : B[idx]-A
		xNew = xB
		opNew = sign ? hcat(B[1:idx-1],S,B[idx+1:N]  ) : -hcat(B[1:idx-1],S,B[idx+1:N]  )
	else
		xNew  = (xA...,xB...)
		opNew = sign ? hcat(A,B) : hcat(A,-B)
	end

	return xNew, opNew
end

#unsigned sum: HCAT+HCAT
function Usum_op(xA::NTuple{NA,Variable},
		 xB::NTuple{NB,Variable},
		 A::L1,
		 B::L2,sign::Bool) where {NA,NB,M,
					  L1<:HCAT{M,NB},
					  L2<:HCAT{M,NB}     }
	xNew = xA
	opNew = A
	for i in eachindex(xB)
		xNew, opNew = Usum_op(xNew, (xB[i],), opNew, B[i], sign)
	end
return xNew,opNew
end

#unsigned sum: multivar AbstractOperator + AbstractOperator
function Usum_op(xA::NTuple{N,Variable},
		 xB::Tuple{Variable},
		 A::AbstractOperator,
		 B::AbstractOperator,sign::Bool) where {N}
	if xB[1] in xA
		Z = Zeros(A)       #this will be an HCAT
		xNew, opNew = Usum_op(xA,xB,Z,B,sign)
		opNew += A
	else
		xNew  = (xA...,xB...)
		opNew = sign ? hcat(A,B) : hcat(A,-B)
	end
	return xNew, opNew
end

"""
`+(ex::AbstractExpression, b::Union{AbstractArray,Number})`

Add a scalar or an `Array` to an expression: 

# Example

```julia
julia> x = Variable(10)
Variable(Float64, (10,))

julia> ex = x+4

```

Notice that in order to add an array to `ex`, `b` must belong to  
of the associate `AbstractOperator` of `ex`. 

```julia
julia> b = randn(10);

julia> size(b), eltype(b)
((10,), Float64)

julia> size(operator(ex),1), codomainType(operator(ex))
((10,), Float64)

julia> ex + b

```

"""
function (+)(a::AbstractExpression, b::Union{AbstractArray,Number}) 
	A = convert(Expression,a)
	return Expression{length(A.x)}(A.x,operator(A),displacement(A)+b)
end

(+)(a::Union{AbstractArray,Number}, b::AbstractExpression) = b+a

function (-)(a::AbstractExpression, b::Union{AbstractArray,Number})
	A = convert(Expression,a)
	return Expression{length(A.x)}(A.x,operator(A),displacement(A)-b)
end

function (-)(a::Union{AbstractArray,Number}, b::AbstractExpression)
	B = convert(Expression,b)
	return Expression{length(B.x)}(B.x,-operator(B),a-displacement(B))
end
# sum with array/scalar

#broadcasted + -
import Base: broadcast

function broadcast(::typeof(+),a::AbstractExpression, b::AbstractExpression)
	A = convert(Expression,a)
	B = convert(Expression,b)
	if size(operator(A),1) != size(operator(B),1)
		if prod(size(operator(A),1)) > prod(size(operator(B),1))
			da = A.d
			db = B.d
			A = Expression{length(A.x)}(A.x,A.L,0.) #remove displacement
			B = Expression{length(B.x)}(variables(B),
						    BroadCast(operator(B),size(operator(A),1)), 
						    0.)
		elseif prod(size(operator(B),1)) > prod(size(operator(A),1))
			da = A.d
			db = B.d
			A = Expression{length(A.x)}(variables(A),
						    BroadCast(operator(A),size(operator(B),1)), 
						    0.)
			B = Expression{length(B.x)}(B.x,B.L,0.) #remove displacement
		end
		return A+B+(da.+db)
	end
	return A+B
end

function broadcast(::typeof(-),a::AbstractExpression, b::AbstractExpression)
	A = convert(Expression,a)
	B = convert(Expression,b)
	if size(operator(A),1) != size(operator(B),1)
		if prod(size(operator(A),1)) > prod(size(operator(B),1))
			da = A.d
			db = B.d
			A = Expression{length(A.x)}(A.x,A.L,0.) #remove displacement
			B = Expression{length(B.x)}(variables(B),
						    BroadCast(operator(B),size(operator(A),1)), 
						    0.)
		elseif prod(size(operator(B),1)) > prod(size(operator(A),1))
			da = A.d
			db = B.d
			A = Expression{length(A.x)}(variables(A),
						    BroadCast(operator(A),size(operator(B),1)), 
						    0.)
			B = Expression{length(B.x)}(B.x,B.L,0.) #remove displacement
		end
		return A-B+(da.-db)
	end
	return A-B
end
