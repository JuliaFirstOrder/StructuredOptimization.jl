
export MyOperator

immutable MyOperator <: LinearOperator
	codomainType::Type
	dim_out::Tuple
	domainType::Type
	dim_in::Tuple
	Fwd!::Function
	Adj!::Function
end

size(L::MyOperator) = (L.dim_out, L.dim_in)

# Constructors
MyOperator(codomainType::Type, dim_out::Int64, 
	   domainType::Type, dim_in::Int64, 
	   Fwd!::Function, Adj!::Function ) =
MyOperator(codomainType, (dim_out,), domainType, (dim_in,), Fwd!, Adj! )

MyOperator(domainType::Type, dim_out::Int64, dim_in::Int64, 
	   Fwd!::Function, Adj!::Function ) =
MyOperator(domainType, (dim_out,), domainType, (dim_in,), Fwd!, Adj! )

MyOperator(domainType::Type, dim_out::Tuple, dim_in::Tuple, 
	   Fwd!::Function, Adj!::Function ) =
MyOperator(domainType, dim_out, domainType, dim_in, Fwd!, Adj! )

MyOperator(dim_out::Int64, dim_in::Int64, 
	   Fwd!::Function, Adj!::Function ) =
MyOperator(Float64, (dim_out,), Float64, (dim_in,), Fwd!, Adj! )

MyOperator(dim_out::Tuple, dim_in::Tuple, 
	   Fwd!::Function, Adj!::Function ) =
MyOperator(Float64, dim_out, Float64, dim_in, Fwd!, Adj! )

MyOperator(dim_out::Tuple, Fwd!::Function, Adj!::Function ) =
MyOperator(Float64, dim_out, Float64, dim_out, Fwd!, Adj! )

MyOperator(dim_out::Int, Fwd!::Function, Adj!::Function ) =
MyOperator(Float64, (dim_out,), Float64, (dim_out,), Fwd!, Adj! )

# Operators
A_mul_B!( y::AbstractArray, L::MyOperator, b::AbstractArray) = L.Fwd!(y,b)
Ac_mul_B!(y::AbstractArray, L::MyOperator, b::AbstractArray) = L.Adj!(y,b)

#Properties
fun_name(L::MyOperator)  = "User Defined Operator"

#utils
