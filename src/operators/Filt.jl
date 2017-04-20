export Filt

immutable Filt <: LinearOperator
	domainType::Type
	dim_in::Tuple
	b::AbstractVector
	a::AbstractVector
	si::AbstractVector
	function Filt(domainType,dim_in,b,a)
		isempty(b) && throw(ArgumentError("filter vector b must be non-empty"))
		isempty(a) && throw(ArgumentError("filter vector a must be non-empty"))
		a[1] == 0  && throw(ArgumentError("filter vector a[1] must be nonzero"))
		as = length(a)
		bs = length(b)
		sz = max(as, bs)
		silen = sz - 1

		# Filter coefficient normalization
		if a[1] != 1
			norml = a[1]
			a ./= norml
			b ./= norml
		end
		# Pad the coefficients with zeros if needed
		bs<sz   && (b = copy!(zeros(eltype(b), sz), b))
		1<as<sz && (a = copy!(zeros(eltype(a), sz), a))

		si = zeros(promote_type(eltype(b), eltype(a)), max(length(a), length(b))-1)
		new(domainType,dim_in,b,a,si)
	end
end

size(L::Filt) = L.dim_in, L.dim_in

# Constructors

Filt{D1}(dim_in::Int,  b::AbstractVector{D1}, a::AbstractVector{D1}) = 
Filt(eltype(b),(dim_in,), b, a)

Filt{D1}(dim_in::Tuple,  b::AbstractVector{D1}, a::AbstractVector{D1}) = 
Filt(eltype(b), dim_in, b, a)

Filt{D1}(x::AbstractArray{D1}, b::AbstractVector{D1}, a::AbstractVector{D1}) = 
Filt(eltype(x),size(x), h)

# Operators

function A_mul_B!{T}(y::AbstractArray{T},L::Filt,x::AbstractArray{T})
	for col = 1:size(x,2)
		length(L.a) != 1 ? iir!(y,L.b,L.a,x,L.si,col,col) : fir!(y,L.b,x,L.si,col,col)
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T},L::Filt,x::AbstractArray{T})
	for col = 1:size(x,2)
		length(L.a) != 1 ? iir_rev!(y,L.b,L.a,x,L.si,col,col) : fir_rev!(y,L.b,x,L.si,col,col)
	end
end

fun_name(L::Filt)  = size(L.a,1) != 1 ? "IIR" : "FIR"

function iir!(y, b, a, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi - a[j+1]*val
        end
        si[silen] = b[silen+1]*xi - a[silen+1]*val
        y[i,coly] = val
    end
    si .= 0. #reset state
end

function iir_rev!(y, b, a, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=size(x, 1):-1:1
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi - a[j+1]*val
        end
        si[silen] = b[silen+1]*xi - a[silen+1]*val
        y[i,coly] = val
    end
    si .= 0.
end

function fir!(y, b, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi
        end
        si[silen] = b[silen+1]*xi
        y[i,coly] = val
    end
    si .= 0.
end

function fir_rev!(y, b, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=size(x, 1):-1:1
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi
        end
        si[silen] = b[silen+1]*xi
        y[i,coly] = val
    end
    si .= 0.
end
