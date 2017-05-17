export FiniteDiff

immutable FiniteDiff{N,D} <: LinearOperator
	domainType::Type
	dim_in::NTuple{N,Int}
end

# Constructors

FiniteDiff(domainType::Type, dim_in::Tuple, dir::Int64 = 1) =
FiniteDiff{length(dim_in),dir}(domainType, dim_in)

FiniteDiff(dim_in::Tuple, dir::Int64 = 1) =
FiniteDiff{length(dim_in),dir}(Float64, dim_in)

FiniteDiff{T,N}(x::AbstractArray{T,N}, dir::Int64 = 1)  = FiniteDiff{N,dir}(eltype(x), size(x))

# Mappings

function A_mul_B!{T}(y::AbstractArray{T,1},L::FiniteDiff{1,1},b::AbstractArray{T,1})
	for l = 1:length(b)
		y[l] = l == 1 ? b[l+1]-b[l] : b[l]-b[l-1]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,1},L::FiniteDiff{1,1},b::AbstractArray{T,1})
	for l = 1:length(b)
		y[l] =
		l == 1 ? -(b[l] + b[l+1]) :
		l == 2 ?   b[l] + b[l-1] - b[l+1] :
		l == length(b) ? b[l] : b[l]-b[l+1]

	end
end

function A_mul_B!{T}(y::AbstractArray{T,2}, L::FiniteDiff{2,1}, b::AbstractArray{T,2})
	for l = 1:size(b,1), m = 1:size(b,2)
		y[l,m] = l == 1 ? b[l+1,m]-b[l,m] : b[l,m]-b[l-1,m]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,2}, L::FiniteDiff{2,1}, b::AbstractArray{T,2})
	for l = 1:size(b,1), m = 1:size(b,2)
		y[l,m] =
		l == 1 ? -(b[l,m] + b[l+1,m]) :
		l == 2 ?   b[l,m] + b[l-1,m] - b[l+1,m] :
		l == size(b,1) ? b[l,m] : b[l,m]-b[l+1,m]
	end
end

function A_mul_B!{T}(y::AbstractArray{T,2},L::FiniteDiff{2,2},b::AbstractArray{T,2})
	for l = 1:size(b,1), m = 1:size(b,2)
		y[l,m] = m == 1 ? b[l,m+1]-b[l,m] : b[l,m]-b[l,m-1]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,2},L::FiniteDiff{2,2},b::AbstractArray{T,2})
	for l = 1:size(b,1), m = 1:size(b,2)
		y[l,m] =
		m == 1 ? -(b[l,m] + b[l,m+1]) :
		m == 2 ?   b[l,m] + b[l,m-1] - b[l,m+1] :
		m == size(b,2) ? b[l,m] : b[l,m]-b[l,m+1]
	end
end

function A_mul_B!{T}(y::AbstractArray{T,3},L::FiniteDiff{3,1},b::AbstractArray{T,3})
	for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3)
		y[l,m,n] = l == 1 ? b[l+1,m,n]-b[l,m,n] : b[l,m,n]-b[l-1,m,n]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,3},L::FiniteDiff{3,1},b::AbstractArray{T,3})
	for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3)
		y[l,m,n] =
		l == 1 ? -(b[l,m,n] + b[l+1,m,n]) :
		l == 2 ?   b[l,m,n] + b[l-1,m,n] - b[l+1,m,n] :
		l == size(b,1) ? b[l,m,n] : b[l,m,n]-b[l+1,m,n]
	end
end

function A_mul_B!{T}(y::AbstractArray{T,3},L::FiniteDiff{3,2},b::AbstractArray{T,3})
	for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3)
		y[l,m,n] = m == 1 ? b[l,m+1,n]-b[l,m,n] : b[l,m,n]-b[l,m-1,n]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,3},L::FiniteDiff{3,2},b::AbstractArray{T,3})
	for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3)
		y[l,m,n] =
		m == 1 ? -(b[l,m,n] + b[l,m+1,n]) :
		m == 2 ?   b[l,m,n] + b[l,m-1,n] - b[l,m+1,n] :
		m == size(b,2) ? b[l,m,n] : b[l,m,n]-b[l,m+1,n]
	end
end

function A_mul_B!{T}(y::AbstractArray{T,3},L::FiniteDiff{3,3},b::AbstractArray{T,3})
	for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3)
		y[l,m,n] = n == 1 ? b[l,m,n+1]-b[l,m,n] : b[l,m,n]-b[l,m,n-1]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,3},L::FiniteDiff{3,3},b::AbstractArray{T,3})
	for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3)
		y[l,m,n] =
		n == 1 ? -(b[l,m,n] + b[l,m,n+1]) :
		n == 2 ?   b[l,m,n] + b[l,m,n-1] - b[l,m,n+1] :
		n == size(b,3) ? b[l,m,n] : b[l,m,n]-b[l,m,n+1]
	end
end

# Properties

size(L::FiniteDiff) = (L.dim_in, L.dim_in)

fun_name{N,D}(L::FiniteDiff{N,D})  = "Finite Diff. Op. dir = $D"
