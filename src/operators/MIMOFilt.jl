export MIMOFilt

immutable MIMOFilt <: LinearOperator
	domainType::Type
	dim_out::Tuple{Int,Int}
	dim_in::Tuple{Int,Int}
	B::Vector{AbstractVector}
	A::Vector{AbstractVector}
	SI::Vector{AbstractVector}

	function MIMOFilt(domainType,dim_in,b,a)

		length(b) != length(a) && throw(ArgumentError("filter vectors b must be as many as a"))

		mod(length(b),dim_in[2]) !=0 && throw(ArgumentError("wrong number of filters"))
		dim_out = (dim_in[1], div(length(b),dim_in[2]) )

		B,A,SI = [],[],[]

		for i = 1:length(b)
			isempty(b[i]) && throw(ArgumentError("filter vector b[$i] must be non-empty"))
			isempty(a[i]) && throw(ArgumentError("filter vector a[$i] must be non-empty"))
			a[i][1] == 0  && throw(ArgumentError("filter vector a[$i][1] must be nonzero"))

			as = length(a[i])
			bs = length(b[i])
			sz = max(as, bs)
			silen = sz - 1


			# Filter coefficient normalization
			if a[i][1] != 1
				norml = a[i][1]
				a[i] ./= norml
				b[i] ./= norml
			end

			# Pad the coefficients with zeros if needed
			bs<sz   && (b[i] = copy!(zeros(eltype(b[i]), sz), b[i]))
			1<as<sz && (a[i] = copy!(zeros(eltype(a[i]), sz), a[i]))

			si = zeros(promote_type(eltype(b[i]),
			   	   eltype(a[i])), max(length(a[i]), length(b[i]))-1)

			push!(B,b[i])
			push!(A,a[i])
			push!(SI,si)
		end
		new(domainType, dim_out, dim_in, B, A, SI)
	end
end

# Constructors

MIMOFilt{D1<:AbstractVector}(dim_in::Tuple,  b::Vector{D1}, a::Vector{D1}) =
MIMOFilt(eltype(b[1]), dim_in, b, a)

MIMOFilt{D1<:AbstractVector}(dim_in::Tuple,  b::Vector{D1}) =
MIMOFilt(eltype(b[1]), dim_in, b, [[1.0] for i in eachindex(b)])

MIMOFilt{D1<:AbstractVector}(x::AbstractMatrix,  b::Vector{D1}, a::Vector{D1}) =
MIMOFilt(eltype(x), size(x), b, a)

MIMOFilt{D1<:AbstractVector}(x::AbstractMatrix,  b::Vector{D1}) =
MIMOFilt(eltype(x), size(x), b, [[1.0] for i in eachindex(b)])

# Mappings

function A_mul_B!{T}(y::AbstractArray{T},L::MIMOFilt,x::AbstractArray{T})
	cnt = 0
	cx  = 0
	y .= 0. #TODO avoid this?
	for cy = 1:L.dim_out[2]
		cnt += 1
		cx  += 1
		length(L.A[cnt]) != 1 ? add_iir!(y,L.B[cnt],L.A[cnt],x,L.SI[cnt],cy,cx) :
		add_fir!(y,L.B[cnt],x,L.SI[cnt],cy,cx)

		for c2 = 2:L.dim_in[2]
			cnt += 1
			cx  += 1
			length(L.A[cnt]) != 1 ? add_iir!(y,L.B[cnt],L.A[cnt],x,L.SI[cnt],cy,cx) :
			add_fir!(y,L.B[cnt],x,L.SI[cnt],cy,cx)
		end
		cx = 0
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T},L::MIMOFilt,x::AbstractArray{T})
	cnt = 0
	cx  = 0
	y .= 0. #TODO avoid this?
	for cy = 1:L.dim_out[2]
		cnt += 1
		cx  += 1
		length(L.A[cnt]) != 1 ? add_iir_rev!(y,L.B[cnt],L.A[cnt],x,L.SI[cnt],cx,cy) :
		add_fir_rev!(y,L.B[cnt],x,L.SI[cnt],cx,cy)

		for c2 = 2:L.dim_in[2]
			cnt += 1
			cx  += 1
			length(L.A[cnt]) != 1 ? add_iir_rev!(y,L.B[cnt],L.A[cnt],x,L.SI[cnt],cx,cy) :
			add_fir_rev!(y,L.B[cnt],x,L.SI[cnt],cx,cy)
		end
		cx = 0
	end
end

# Properties

size(L::MIMOFilt) = L.dim_out, L.dim_in

fun_name(L::MIMOFilt)  = "MIMO Filt"

# Utilities

function add_iir!(y, b, a, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi - a[j+1]*val
        end
        si[silen] = b[silen+1]*xi - a[silen+1]*val
        y[i,coly] += val
    end
    si .= 0. #reset state
end

function add_iir_rev!(y, b, a, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=size(x, 1):-1:1
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi - a[j+1]*val
        end
        si[silen] = b[silen+1]*xi - a[silen+1]*val
        y[i,coly] += val
    end
    si .= 0.
end

function add_fir!(y, b, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi
        end
        si[silen] = b[silen+1]*xi
        y[i,coly] += val
    end
    si .= 0.
end

function add_fir_rev!(y, b, x, si, coly, colx)
    silen = length(si)
    @inbounds for i=size(x, 1):-1:1
        xi = x[i,colx]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi
        end
        si[silen] = b[silen+1]*xi
        y[i,coly] += val
    end
    si .= 0.
end
