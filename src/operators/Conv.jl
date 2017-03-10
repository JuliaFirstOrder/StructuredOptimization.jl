import Base: conv, xcorr

immutable Conv{D1,D2} <: LinearOperator{D1,D2}
	x::OptVar
	h::AbstractVector
	dim::Tuple
	isTranspose::Bool
end
size(A::Conv) = A.dim 

conv{D1}(x::OptVar{D1},h::AbstractVector) = 
Conv{D1,D1}(x, h, (size(x), ( length(x.x)+length(h)-1, )  ), false)
conv{D1}(h::AbstractVector,x::OptVar{D1}) = 
Conv{D1,D1}(x, h, (size(x), ( length(x.x)+length(h)-1, )  ), false)

function A_mul_B!{T}(y::AbstractVector{T},A::Conv,b::AbstractVector{T})
	if A.isTranspose
		l =floor(Int64,size(A,1)[1]/2)
		idx = l+1:l+size(A,2)[1]
		y .= conv(b,A.h)[idx]
	else
		y .= conv(A.h,b)
	end
end

transpose{D1}( A::Conv{D1,D1} ) = Xcorr{D1,D1}(A.x, A.h, (A.dim[2],A.dim[1]), true) 

fun_name(A::Conv)  = "Convolution Operator"

conv(B::LinearOperator, args...) = NestedLinearOperator(conv, B, args...)


immutable Xcorr{D1,D2} <: LinearOperator{D1,D2}
	x::OptVar
	h::AbstractVector
	dim::Tuple
	isTranspose::Bool
end
size(A::Xcorr) = A.dim

xcorr{D1}(x::OptVar{D1}, h::AbstractVector) = 
Xcorr{D1,D1}(x,h,(size(x), ( 2*max(length(x.x),length(h))-1, )  ), false)

function A_mul_B!{T}(y::AbstractVector{T},A::Xcorr,b::AbstractVector{T})
	if A.isTranspose
		y .= xcorr(b,A.h)[size(A,1)[1]:end-length(A.h)+1]
	else
		y .= xcorr(A.h,b)
	end
end

transpose{D1}(A::Xcorr{D1,D1} ) = Conv{D1,D1}(A.x,A.h,(A.dim[2],A.dim[1]), true) 

fun_name(A::Xcorr)  = "Cross Correlation Operator"

xcorr(B::LinearOperator, args...) = NestedLinearOperator(xcorr, B, args...)




