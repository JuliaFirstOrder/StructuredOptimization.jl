import Base: conv, xcorr
export Conv, Xcorr

immutable Conv{D1,D2} <: LinearOperator{D1,D2}
	h::AbstractVector
	dim::Tuple
	isTranspose::Bool
end
size(A::Conv) = A.dim 

Conv{D1}(h::AbstractVector{D1}, x::AbstractVector{D1}) = 
Conv{D1,D1}(h, (size(x), ( length(x)+length(h)-1, )  ), false)

function conv(h::AbstractVector, x::OptVar) 
	A = Conv(h,x.x) 
	Affine([x], A, A', Nullable{Vector{AbstractArray}}() )
end

conv(x::OptVar,h::AbstractVector) = conv(h,x)

function A_mul_B!{T}(y::AbstractVector{T},A::Conv,b::AbstractVector{T})
	if A.isTranspose
		l =floor(Int64,size(A,1)[1]/2)
		idx = l+1:l+size(A,2)[1]
		y .= conv(b,A.h)[idx]
	else
		y .= conv(A.h,b)
	end
end

transpose{D1}( A::Conv{D1,D1} ) = Xcorr{D1,D1}(A.h, (A.dim[2],A.dim[1]), true) 

fun_name(A::Conv)  = "Convolution Operator"

conv(B::AffineOperator, args...) = NestedLinearOperator(conv, B, args...)

immutable Xcorr{D1,D2} <: LinearOperator{D1,D2}
	h::AbstractVector
	dim::Tuple
	isTranspose::Bool
end
size(A::Xcorr) = A.dim

Xcorr{D1}(x::AbstractVector{D1}, h::AbstractVector{D1}) = 
Xcorr{D1,D1}(h,(size(x), ( 2*max(length(x),length(h))-1, )  ), false)

function xcorr(x::OptVar, h::AbstractVector) 
	A = Xcorr(x.x,h) 
	Affine([x], A, A', Nullable{Vector{AbstractArray}}() )
end

function A_mul_B!{T}(y::AbstractVector{T},A::Xcorr,b::AbstractVector{T})
	if A.isTranspose
		y .= xcorr(b,A.h)[size(A,1)[1]:end-length(A.h)+1]
	else
		y .= xcorr(A.h,b)
	end
end

transpose{D1}(A::Xcorr{D1,D1} ) = Conv{D1,D1}(A.h,(A.dim[2],A.dim[1]), true) 

fun_name(A::Xcorr)  = "Cross Correlation Operator"

xcorr(B::AffineOperator, args...) = NestedLinearOperator(xcorr, B, args...)




