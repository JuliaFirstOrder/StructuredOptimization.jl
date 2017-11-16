import Base: +, -, *, convert
export variables, operator, displacement

immutable Expression{N} <: AbstractExpression
	x::NTuple{N,Variable}
	L::AbstractOperator
	d::Union{Number, AbstractArray}
	function Expression{N}(x::NTuple{N,Variable},L,d) where {N}

		# checks on L
		ndoms(L,1) > 1 && throw(ArgumentError(
	"cannot create expression with LinearOperator with `ndoms(L,1) > 1`"))

		#checks on x
		szL = size(L,2)
		szx = size.(x)
		check_sz = length(szx) == 1 ? szx[1] != szL : szx != szL
		check_sz && throw(ArgumentError(
	"Size of the operator domain $(size(L, 2)) must match size of the variable $(size.(x))"))

		dmL = domainType(L)
		dmx = eltype.(x)
		check_dm = length(dmx) == 1 ? dmx[1] != dmL : dmx != dmL
		check_dm && throw(ArgumentError(
	"Type of the operator domain $(domainType(L)) must match type of the variable $(eltype.(x))"))

		#checks on d
		if typeof(d) <: AbstractArray
			size(L,1) != size(d) && throw(ArgumentError(
	"cannot sum Array of dimensions $(size(d)) with expression with codomain size $(size(L,1))"))
			codomainType(L) != eltype(d) && throw(ArgumentError(
	"cannot sum $(typeof(d)) with expression with codomain $(codomainType(L))"))
		else
			if typeof(d) <: Complex && codomainType(L) <: Real
				throw(ArgumentError(
	"cannot sum $(typeof(d)) to an expression with codomain $(codomainType(L))"))
			end
		end
		new{N}(x,L,d)
	end
end

convert{T,N,A}(::Type{Expression},x::Variable{T,N,A}) =
Expression{1}((x,),Eye(T,size(x)),zero(T))

# constructors

# multipy expressions with AbstractOperator
function (*)(L::T1, a::T2) where {T1 <: AbstractOperator, T2 <: AbstractExpression}
	A = convert(Expression,a)
	if typeof(displacement(A)) <: Number
		d = displacement(A) == 0. ? zero(codomainType(L)) :
		L*(displacement(A)*ones(codomainType(operator(A)),size(operator(A),1)))
	else
		d = L*displacement(A)
	end
	Expression{length(A.x)}(A.x,L*operator(A),d)
end

# sum expressions
function (+)(a::T1, b::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression}
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

function (-)(a::T1, b::T2) where {T1 <: AbstractExpression, T2 <: AbstractExpression}
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
		 A::L1,
		 B::L2,sign::Bool) where {L1<:AbstractOperator,
					  L2<:AbstractOperator}
	xNew  = (xA...,xB...)
	opNew = sign ? hcat(A,B) : hcat(A,-B)
	return xNew, opNew
end

#unsigned sum: HCAT + AbstractOperator
function Usum_op(xA::NTuple{N,Variable},
		 xB::Tuple{Variable},
		 A::L1,
		 B::L2,sign::Bool) where {N,M,L1<:HCAT{M,N},
					  L2<:AbstractOperator}
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
		 A::L1,
		 B::L2,sign::Bool) where {N,M,
					  L1<:AbstractOperator,
					  L2<:HCAT{M,N}     }
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
		 A::L1,
		 B::L2,sign::Bool) where {N,
					  L1<:AbstractOperator,
					  L2<:AbstractOperator}
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

# sum with array/scalar
function (+)(a::T1, b::T2) where {T1 <: AbstractExpression, T2 <: Union{AbstractArray,Number}}
	A = convert(Expression,a)
	return Expression{length(A.x)}(A.x,operator(A),displacement(A)+b)
end

(+)(a::T1, b::T2) where {T1 <: Union{AbstractArray,Number}, T2 <: AbstractExpression} = b+a

function (-)(a::T1, b::T2) where {T1 <: AbstractExpression, T2 <: Union{AbstractArray,Number}}
	A = convert(Expression,a)
	return Expression{length(A.x)}(A.x,operator(A),displacement(A)-b)
end

function (-)(a::T1, b::T2) where {T1 <: Union{AbstractArray,Number}, T2 <: AbstractExpression}
	B = convert(Expression,b)
	return Expression{length(B.x)}(B.x,-operator(B),a-displacement(B))
end

#broadcasted + -
import Base: broadcast

function broadcast(::typeof(+),a::T1,b::T2) where  {T1 <: AbstractExpression, T2 <: AbstractExpression}
	A = convert(Expression,a)
	B = convert(Expression,b)
	if size(operator(A),1) != size(operator(B),1)
		if ndims(operator(A),1) > ndims(operator(B),1) || size(operator(B),1) == (1,)
			da = A.d
			db = B.d
			A = Expression{length(A.x)}(A.x,A.L,0.) #remove displacement
			B = Expression{length(B.x)}(variables(B),
						    BroadCast(operator(B),size(operator(A),1)), 
						    0.)
		elseif ndims(operator(B),1) > ndims(operator(A),1) || size(operator(A),1) == (1,)
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

function broadcast(::typeof(-),a::T1,b::T2) where  {T1 <: AbstractExpression, T2 <: AbstractExpression}
	A = convert(Expression,a)
	B = convert(Expression,b)
	if size(operator(A),1) != size(operator(B),1)
		if ndims(operator(A),1) > ndims(operator(B),1) || size(operator(B),1) == (1,)
			da = A.d
			db = B.d
			A = Expression{length(A.x)}(A.x,A.L,0.) #remove displacement
			B = Expression{length(B.x)}(variables(B),
						    BroadCast(operator(B),size(operator(A),1)), 
						    0.)
		elseif ndims(operator(B),1) > ndims(operator(A),1) || size(operator(A),1) == (1,)
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

# AbstractOperators binding
# special cases
import Base: *, reshape

#NonLinearCompose
function (*)(ex1::T1, ex2::T2) where {T1<:AbstractExpression,T2<:AbstractExpression}
	ex1 = convert(Expression,ex1)
	ex2 = convert(Expression,ex2)
	op = NonLinearCompose(operator(ex1),operator(ex2))
	d = (displacement(ex1) != 0 && displacement(ex2) != 0) ? displacement(ex1)*displacement(ex2) : 
	zero(codomainType(op))
	x = (variables(ex1)...,variables(ex2)...) 
	exp3 = Expression{length(x)}(x,op,d)
	if displacement(ex2) != 0.
		exp3 += Expression{length(ex1.x)}(ex1.x,ex1.L,0.)*displacement(ex2)
	end
	if displacement(ex1) != 0.
		exp3 += displacement(ex1)*Expression{length(ex2.x)}(ex2.x,ex2.L,0.)
	end
	return exp3 
end

#LMatrixOp
function (*)(m::T, a::Union{AbstractVector,AbstractMatrix}) where {T<:AbstractExpression}
	M = convert(Expression,m)
	op = LMatrixOp(codomainType(operator(M)),size(operator(M),1),a)
	return op*M
end

#MatrixOp
function (*)(M::AbstractMatrix, a::T) where {T<:AbstractExpression}
	A = convert(Expression,a)
	op = MatrixOp(codomainType(operator(A)),size(operator(A),1),M)
	return op*A
end

#DiagOp
function Base.broadcast(::typeof(*), d::AbstractArray, a::T) where {T<:AbstractExpression}
	A = convert(Expression,a)
	op = DiagOp(codomainType(operator(A)),size(operator(A),1),d)
	return op*A
end

#Scale
function (*)(coeff::T1, a::T) where {T1<:Number, T<:AbstractExpression}
	A = convert(Expression,a)
	return Expression{length(A.x)}(A.x,coeff*operator(A),coeff*displacement(A))
end

#Reshape
function reshape(a::T, dims...) where {T<:AbstractExpression}
	A = convert(Expression,a)
	op = Reshape(A.L, dims...)
	if typeof(displacement(A)) <: Number
		d = displacement(A)
	else
		d = reshape(displacement(A), dims...)
	end
	return Expression{length(A.x)}(A.x,op,d)
end

imported = [:getindex :GetIndex;
	    :fft      :DFT;
	    :rfft     :RDFT;
	    :irfft    :IRDFT;
	    :ifft     :IDFT;
	    :dct      :DCT;
	    :idct     :IDCT;
	    :conv     :Conv;
	    :xcorr    :Xcorr;
	    :filt     :Filt;
	    ]

exported = [:finitediff :FiniteDiff;
	    :variation  :Variation;
	    :mimofilt   :MIMOFilt;
	    :zeropad    :ZeroPad;
	    :sigmoid    :Sigmoid;
	    :σ          :Sigmoid; #alias
	    ]

#importing functions from Base
for f in  imported[:,1]
	@eval begin
		import Base: $f
	end
end
#exporting functions
for f in  exported[:,1]
	@eval begin
		export $f
	end
end


fun = [imported; exported]
for i = 1:size(fun,1)
	f,fAbsOp = fun[i,1],fun[i,2]
	@eval begin
		function $f{T<:AbstractExpression}(a::T,args...)
			A = convert(Expression,a)
			op = $fAbsOp(codomainType(operator(A)),size(operator(A),1), args...)
			return op*A
		end
	end
end

# Properties
variables(A::Expression)    = A.x
operator(A::Expression)     = A.L
displacement(A::Expression) = A.d
