export FiniteDiff, finiteDiff

immutable FiniteDiff{D1,D2,N,D} <: LinearOperator{D1,D2}
	sign::Bool
	isTranspose::Bool
	dim::Tuple

	FiniteDiff(      isTranspose, dim) = FiniteDiff{D1,D2,N,D}(true, isTranspose, dim)
	FiniteDiff(sign, isTranspose, dim) = new(sign, isTranspose, dim)
end

size(A::FiniteDiff) = (A.dim, A.dim)
-{D1,D2,N,D}(A::FiniteDiff{D1,D2,N,D}) = 
FiniteDiff{D1,D2,N,D}(false == sign(A), A.isTranspose, A.dim) 

FiniteDiff{T,N}(x::AbstractArray{T,N}, dir::Int64 = 1)  = FiniteDiff{T,T,N,dir}(false, size(x))

finiteDiff(x::Variable, dir::Int64 = 1) = (FiniteDiff(~x, dir))*x

transpose{D1,N,D}(A::FiniteDiff{D1,D1,N,D}) = 
FiniteDiff{D1,D1,N,D}(A.sign, !(A.isTranspose), A.dim)

fun_name{D1,D2,N,D}(A::FiniteDiff{D1,D2,N,D})  = "Finite Diff. Op. dir = $D"

function uA_mul_B!{T}(y::AbstractArray{T,1},A::FiniteDiff{T,T,1,1},b::AbstractArray{T,1})
	if A.isTranspose == false 
		for l = 1:length(b)
			y[l] = l == 1 ? b[l+1]-b[l] : b[l]-b[l-1]
		end
	else
		for l = 1:length(b) 
			y[l] = 
			 l == 1 ? -(b[l] + b[l+1]) :
			 l == 2 ?   b[l] + b[l-1] - b[l+1] : 
			 l == length(b) ? b[l] : b[l]-b[l+1]
		end

	end
end

function uA_mul_B!{T}(y::AbstractArray{T,2},A::FiniteDiff{T,T,2,1},b::AbstractArray{T,2})
	if A.isTranspose == false 
		for l = 1:size(b,1), m = 1:size(b,2) 
			y[l,m] = l == 1 ? b[l+1,m]-b[l,m] : b[l,m]-b[l-1,m]
		end
	else
		for l = 1:size(b,1), m = 1:size(b,2) 
			y[l,m] = 
			 l == 1 ? -(b[l,m] + b[l+1,m]) :
			 l == 2 ?   b[l,m] + b[l-1,m] - b[l+1,m] : 
			 l == size(b,1) ? b[l,m] : b[l,m]-b[l+1,m]
		end

	end
end

function uA_mul_B!{T}(y::AbstractArray{T,2},A::FiniteDiff{T,T,2,2},b::AbstractArray{T,2})
	if A.isTranspose == false 
		for l = 1:size(b,1), m = 1:size(b,2) 
			y[l,m] = m == 1 ? b[l,m+1]-b[l,m] : b[l,m]-b[l,m-1]
		end
	else
		for l = 1:size(b,1), m = 1:size(b,2) 
			y[l,m] = 
			 m == 1 ? -(b[l,m] + b[l,m+1]) :
			 m == 2 ?   b[l,m] + b[l,m-1] - b[l,m+1] : 
			 m == size(b,2) ? b[l,m] : b[l,m]-b[l,m+1]
		end

	end
end

function uA_mul_B!{T}(y::AbstractArray{T,3},A::FiniteDiff{T,T,3,1},b::AbstractArray{T,3})
	if A.isTranspose == false 
		for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3) 
			y[l,m,n] = l == 1 ? b[l+1,m,n]-b[l,m,n] : b[l,m,n]-b[l-1,m,n]
		end
	else
		for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3) 
			y[l,m,n] = 
			 l == 1 ? -(b[l,m,n] + b[l+1,m,n]) :
			 l == 2 ?   b[l,m,n] + b[l-1,m,n] - b[l+1,m,n] : 
			 l == size(b,1) ? b[l,m,n] : b[l,m,n]-b[l+1,m,n]
		end

	end
end

function uA_mul_B!{T}(y::AbstractArray{T,3},A::FiniteDiff{T,T,3,2},b::AbstractArray{T,3})
	if A.isTranspose == false 
		for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3) 
			y[l,m,n] = m == 1 ? b[l,m+1,n]-b[l,m,n] : b[l,m,n]-b[l,m-1,n]
		end
	else
		for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3) 
			y[l,m,n] = 
			 m == 1 ? -(b[l,m,n] + b[l,m+1,n]) :
			 m == 2 ?   b[l,m,n] + b[l,m-1,n] - b[l,m+1,n] : 
			 m == size(b,2) ? b[l,m,n] : b[l,m,n]-b[l,m+1,n]
		end

	end
end

function uA_mul_B!{T}(y::AbstractArray{T,3},A::FiniteDiff{T,T,3,3},b::AbstractArray{T,3})
	if A.isTranspose == false 
		for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3) 
			y[l,m,n] = n == 1 ? b[l,m,n+1]-b[l,m,n] : b[l,m,n]-b[l,m,n-1]
		end
	else
		for l = 1:size(b,1), m = 1:size(b,2), n = 1:size(b,3) 
			y[l,m,n] = 
			 n == 1 ? -(b[l,m,n] + b[l,m,n+1]) :
			 n == 2 ?   b[l,m,n] + b[l,m,n-1] - b[l,m,n+1] : 
			 n == size(b,3) ? b[l,m,n] : b[l,m,n]-b[l,m,n+1]
		end

	end
end


