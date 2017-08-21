import Base: +, -, *, convert
export variables, operator, displacement

immutable AffineExpression{N} <: AbstractExpression
	x::NTuple{N,Variable}
	L::AbstractOperator
	b::Union{Number, AbstractArray}
	function AffineExpression{N}(x::NTuple{N,Variable},L,b) where {N}

		# checks operator is Linear
		is_linear(L) == false && throw(ArgumentError(
	"AffineExperssion must be linear"))

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

		#checks on b
		if typeof(b) <: AbstractArray
			size(L,1) != size(b) && throw(ArgumentError(
	"cannot sum Array of dimensions $(size(b)) with expression with codomain size $(size(L,1))"))
			codomainType(L) != eltype(b) && throw(ArgumentError(
	"cannot sum $(typeof(b)) with expression with codomain $(codomainType(L))"))
		else
			if typeof(b) <: Complex && codomainType(L) <: Real
				throw(ArgumentError(
	"cannot sum $(typeof(b)) to an expression with codomain $(codomainType(L))"))
			end
		end
		new{N}(x,L,b)
	end
end

convert{T,N,A}(::Type{AffineExpression},x::Variable{T,N,A}) =
AffineExpression{1}((x,),Eye(T,size(x)),zero(T))

# constructors

# multipy expressions
function (*){T1 <: AbstractOperator, T2 <: AbstractExpression}(L::T1, a::T2)
	A = convert(AffineExpression,a)
	if typeof(displacement(A)) <: Number
		b = displacement(A) == 0. ? zero(codomainType(L)) :
		L*(displacement(A)*ones(codomainType(operator(A)),size(operator(A),1)))
	else
		b = L*displacement(A)
	end
	AffineExpression{length(A.x)}(A.x,L*operator(A),b)
end

# sum expressions
function (+){T1 <: AbstractExpression, T2 <: AbstractExpression}(a::T1, b::T2)
	A = convert(AffineExpression,a)
	B = convert(AffineExpression,b)
	b = displacement(A)+displacement(B)
	if variables(A) == variables(B)
		return AffineExpression{length(A.x)}(A.x,operator(A)+operator(B),b)
	else
		opA = operator(A)
		xA = variables(A)
		opB = operator(B)
		xB = variables(B)

		xNew, opNew = Usum_op(xA,xB,opA,opB,true)
		return AffineExpression{length(xNew)}(xNew,opNew,b)
	end

end

function (-){T1 <: AbstractExpression, T2 <: AbstractExpression}(a::T1, b::T2)
	A = convert(AffineExpression,a)
	B = convert(AffineExpression,b)
	b = displacement(A)-displacement(B)
	if variables(A) == variables(B)
		return AffineExpression{length(A.x)}(A.x,operator(A)-operator(B),b)
	else
		opA = operator(A)
		xA = variables(A)
		opB = operator(B)
		xB = variables(B)

		xNew, opNew = Usum_op(xA,xB,opA,opB,false)
		return AffineExpression{length(xNew)}(xNew,opNew,b)
	end

end

function Usum_op{L1<:AbstractOperator,
		 L2<:AbstractOperator}(xA::Tuple{Variable},
				       xB::Tuple{Variable},
				       A::L1,
				       B::L2,sign::Bool)
	xNew  = (xA...,xB...)
	opNew = sign ? hcat(A,B) : hcat(A,-B)
	return xNew, opNew
end

function Usum_op{N,M,L1<:HCAT{M,N},
		 L2<:AbstractOperator}(xA::NTuple{N,Variable},
				       xB::Tuple{Variable},
				       A::L1,
				       B::L2,sign::Bool)
	if xB[1] in xA
		idx = findfirst(xB.==xA)
		S = sign ? A[idx]+B : A[idx]-B
		xNew = xA
		opNew = hcat(A[1:idx-1],S,A[idx+1:N]  )
	else
		xNew  = (xA...,xB...)
		opNew = sign ? hcat(A,B) : hcat(A,-B)
	end

	return xNew, opNew
end

function Usum_op{N,M,L1<:AbstractOperator,
		 L2<:HCAT{M,N}     }(xA::Tuple{Variable},
		                     xB::NTuple{N,Variable},
				     A::L1,
				     B::L2,sign::Bool)
	if xA[1] in xB
		idx = findfirst(xA.==xB)
		S = sign ? A+B[idx] : B[idx]-A
		xNew = xB
		opNew = sign ? hcat(B[1:idx-1],S,B[idx+1:N]  ) : -hcat(B[1:idx-1],S,B[idx+1:N]  )
	else
		xNew  = (xA...,xB...)
		opNew = sign ? hcat(A,B) : hcat(A,-B)
	end

	return xNew, opNew
end

function Usum_op{NA,NB,M,L1<:HCAT{M,NB},
		 L2<:HCAT{M,NB}     }(xA::NTuple{NA,Variable},
		                      xB::NTuple{NB,Variable},
				      A::L1,
				      B::L2,sign::Bool)
	xNew = xA
	opNew = A
	for i in eachindex(xB)
		xNew, opNew = Usum_op(xNew, (xB[i],), opNew, B[i], sign)
	end
	return xNew,opNew
end

# sum with array/scalar
function (+){T1 <: AbstractExpression, T2 <: Union{AbstractArray,Number}}(a::T1, b::T2)
	A = convert(AffineExpression,a)
	return AffineExpression{length(A.x)}(A.x,operator(A),displacement(A)+b)
end

(+){T1 <: Union{AbstractArray,Number}, T2 <: AbstractExpression}(a::T1, b::T2) = b+a

function (-){T1 <: AbstractExpression, T2 <: Union{AbstractArray,Number}}(a::T1, b::T2)
	A = convert(AffineExpression,a)
	return AffineExpression{length(A.x)}(A.x,operator(A),displacement(A)-b)
end

function (-){T1 <: Union{AbstractArray,Number}, T2 <: AbstractExpression}(a::T1, b::T2)
	B = convert(AffineExpression,b)
	return AffineExpression{length(B.x)}(B.x,-operator(B),a-displacement(B))
end

# AbstractOperators binding
# special cases
import Base: *, reshape

#MatrixOp
function (*){T<:AbstractExpression}(M::AbstractMatrix, a::T)
	A = convert(AffineExpression,a)
	op = MatrixOp(codomainType(operator(A)),size(operator(A),1),M)
	return op*A
end

#DiagOp
function Base.broadcast{T<:AbstractExpression}(::typeof(*), d::AbstractArray, a::T)
	A = convert(AffineExpression,a)
	op = DiagOp(codomainType(operator(A)),size(operator(A),1),d)
	return op*A
end

#Scale
function (*){T1<:Number, T<:AbstractExpression}(coeff::T1, a::T)
	A = convert(AffineExpression,a)
	return AffineExpression{length(A.x)}(A.x,coeff*operator(A),displacement(A)-b)
end

#TODO Reshape

imported = [:getindex :GetIndex;
	    :fft      :DFT;
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
			A = convert(AffineExpression,a)
			op = $fAbsOp(codomainType(operator(A)),size(operator(A),1), args...)
			return op*A
		end
	end
end

# Properties
variables(A::AffineExpression)    = A.x
operator(A::AffineExpression)     = A.L
displacement(A::AffineExpression) = A.b
