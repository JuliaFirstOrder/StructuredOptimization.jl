export MyOperator

immutable MyOperator{N,M,C,D} <: LinearOperator
	dim_out::NTuple{N,Int}
	dim_in::NTuple{M,Int}
	Fwd!::Function
	Adj!::Function
end

# Constructors

MyOperator(codomainType::Type, dim_out::Int64,
	   domainType::Type, dim_in::Int64,
	   Fwd!::Function, Adj!::Function ) =
MyOperator{1,1,codomainType,domainType}((dim_out,), (dim_in,), Fwd!, Adj! )

MyOperator(domainType::Type, dim_out::Int64, dim_in::Int64,
	   Fwd!::Function, Adj!::Function ) =
MyOperator{1,1,domainType,domainType}((dim_out,), (dim_in,), Fwd!, Adj! )

MyOperator(domainType::Type, dim_out::Tuple, dim_in::Tuple,
	   Fwd!::Function, Adj!::Function ) =
MyOperator{length(dim_out),length(dim_in), domainType, domainType}(dim_out, dim_in, Fwd!, Adj! )

MyOperator(dim_out::Int64, dim_in::Int64,
	   Fwd!::Function, Adj!::Function ) =
MyOperator{1,1,Float64,Float64}((dim_out,), (dim_in,), Fwd!, Adj! )

MyOperator(dim_out::Tuple, dim_in::Tuple,
	   Fwd!::Function, Adj!::Function ) =
MyOperator{length(dim_out),length(dim_in),Float64,Float64}(dim_out, dim_in, Fwd!, Adj! )

MyOperator(dim_out::Tuple, Fwd!::Function, Adj!::Function ) =
MyOperator{length(dim_out),length(dim_out),Float64,Float64}(dim_out, dim_out, Fwd!, Adj! )

MyOperator(dim_out::Int, Fwd!::Function, Adj!::Function ) =
MyOperator{1,1,Float64,Float64}((dim_out,), (dim_out,), Fwd!, Adj! )

# Mappings

A_mul_B!{N,M,C,D}( y::Array{C,N}, L::MyOperator{N,M,C,D}, b::Array{D,M}) = L.Fwd!(y,b)
Ac_mul_B!{N,M,C,D}(y::Array{C,N}, L::MyOperator{N,M,C,D}, b::Array{D,M}) = L.Adj!(y,b)

# Properties

size(L::MyOperator) = (L.dim_out, L.dim_in)

codomainType{N,M,C,D}(L::MyOperator{N,M,C,D}) = C
  domainType{N,M,C,D}(L::MyOperator{N,M,C,D}) = D

fun_name(L::MyOperator)  = "User Defined Operator"
