export Conv, Xcorr

immutable Conv{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	h::AbstractVector
	dim::Tuple
	isTranspose::Bool
	Conv(sign,h,dim,isTranspose) = new(sign,h,dim,isTranspose)
	Conv(     h,dim,isTranspose) = new(true,h,dim,isTranspose)
end

size(A::Conv) = A.dim
-{D1,D2}(A::Conv{D1,D2}) = Conv{D1,D2}(false == sign(A), A.h, A.dim, A.isTranspose)

Conv{D1}(h::AbstractVector{D1}, x::AbstractVector{D1}) =
Conv{D1,D1}(h, (size(x), ( length(x)+length(h)-1, )  ), false)

function uA_mul_B!{T}(y::AbstractVector{T},A::Conv,b::AbstractVector{T})
	if A.isTranspose
		l =floor(Int64,size(A,1)[1]/2)
		idx = l+1:l+size(A,2)[1]
		y .= conv(b,A.h)[idx]
	else
		y .= conv(A.h,b)
	end
end

transpose{D1}( A::Conv{D1,D1} ) = Xcorr{D1,D1}(sign(A), A.h, (A.dim[2],A.dim[1]), !(A.isTranspose))

immutable Xcorr{D1,D2} <: LinearOperator{D1,D2}
	sign::Bool
	h::AbstractVector
	dim::Tuple
	isTranspose::Bool

	Xcorr(sign,h,dim,isTranspose) = new(sign,h,dim,isTranspose)
	Xcorr(     h,dim,isTranspose) = new(true,h,dim,isTranspose)
end
size(A::Xcorr) = A.dim
-{D1,D2}(A::Xcorr{D1,D2}) = Xcorr{D1,D2}(false == sign(A), A.h, A.dim, A.isTranspose)

Xcorr{D1}(x::AbstractVector{D1}, h::AbstractVector{D1}) =
Xcorr{D1,D1}(h,(size(x), ( 2*max(length(x),length(h))-1, )  ), false)

function uA_mul_B!{T}(y::AbstractVector{T},A::Xcorr,b::AbstractVector{T})
	if A.isTranspose
		y .= xcorr(b,A.h)[size(A,1)[1]:end-length(A.h)+1]
	else
		y .= xcorr(A.h,b)
	end
end

transpose{D1}(A::Xcorr{D1,D1} ) = Conv{D1,D1}(sign(A),A.h,(A.dim[2],A.dim[1]), !(A.isTranspose))

fun_name(A::Conv)  = "Convolution Operator"
fun_name(A::Xcorr)  = "Cross Correlation Operator"

################################################################################
# FROM HERE ON IT IS USERS' SYNTAX
################################################################################

import Base: conv, xcorr

function conv(h::AbstractVector, x::Variable)
	A = Conv(h,x.x)
	Affine([x], A, A', Nullable{AbstractArray}() )
end

conv(x::Variable,h::AbstractVector) = conv(h,x)

conv(B::AffineOperator, args...) = NestedLinearOperator(conv, B, args...)

function xcorr(x::Variable, h::AbstractVector)
	A = Xcorr(x.x,h)
	Affine([x], A, A', Nullable{AbstractArray}() )
end

xconv(x::Variable,h::AbstractVector) = xconv(h,x)

xcorr(B::AffineOperator, args...) = NestedLinearOperator(xcorr, B, args...)
