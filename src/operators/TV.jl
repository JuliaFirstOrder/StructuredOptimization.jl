export tv, TV

immutable TV{D1,D2,N} <: LinearOperator{D1,D2}
	sign::Bool
	isTranspose::Bool
	dim::Tuple

	TV(      isTranspose, dim) = TV{D1,D2,N}(true, isTranspose, dim)
	TV(sign, isTranspose, dim) = new(sign, isTranspose, dim)
end

size(A::TV) = (A.dim[1], A.dim[2])
-{D1,D2,N}(A::TV{D1,D2,N}) = 
TV{D1,D2,N}(false == sign(A), A.isTranspose, A.dim) 

TV{T,N}(x::AbstractArray{T,N})  = TV{T,T,N}(false, (size(x),(prod(size(x)),ndims(x)) ))

tv(x::Variable) = (TV(~x))*x

transpose{D1,N}(A::TV{D1,D1,N}) = 
TV{D1,D1,N}(A.sign, !(A.isTranspose), (A.dim[2],A.dim[1]))

fun_name(A::TV)  = "Total Variation"

function uA_mul_B!{T}(y::AbstractArray{T,2},A::TV{T,T,2},b::AbstractArray{T,2})
	if A.isTranspose == false 
		cnt = 0
		for m = 1:size(b,2), l = 1:size(b,1) 
			cnt += 1
			y[cnt,1] = l == 1 ? b[l+1,m]-b[l,m] : b[l,m]-b[l-1,m]
			y[cnt,2] = m == 1 ? b[l,m+1]-b[l,m] : b[l,m]-b[l,m-1]
		end
	else
		cnt = 0
		Nx, Ny = size(y,1), size(y,2)
		for m = 1:Ny, l = 1:Nx 
			cnt += 1
			y[l,m] = (( 
			 l == 1  ? -(b[cnt,1] + b[cnt+1,1]) :
			 l == 2  ?   b[cnt,1] + b[cnt-1,1] - b[cnt+1,1] : 
			 l == Nx ?   b[cnt,1] : b[cnt,1]   - b[cnt+1,1])
			 +(
			 m == 1  ? -(b[cnt,2] + b[cnt+Nx,2]) :
			 m == 2  ?   b[cnt,2] + b[cnt-Nx,2] - b[cnt+Nx,2] : 
			 m == Ny ?   b[cnt,2] : b[cnt,2]    - b[cnt+Nx,2] ))
		end

	end
end

function uA_mul_B!{T}(y::AbstractArray,A::TV{T,T,3},b::AbstractArray)
	if A.isTranspose == false 
		cnt = 0
		for n = 1:size(b,3), m = 1:size(b,2), l = 1:size(b,1) 
			cnt += 1
			y[cnt,1] = l == 1 ? b[l+1,m,n]-b[l,m,n] : b[l,m,n]-b[l-1,m,n]
			y[cnt,2] = m == 1 ? b[l,m+1,n]-b[l,m,n] : b[l,m,n]-b[l,m-1,n]
			y[cnt,3] = n == 1 ? b[l,m,n+1]-b[l,m,n] : b[l,m,n]-b[l,m,n-1]
		end
	else
		cnt = 0
		Nx, Ny, Nz = size(y,1), size(y,2), size(y,3)
		Nxy = Nx*Ny
		for n = 1:Nz, m = 1:Ny, l = 1:Nx 
			cnt += 1
			y[l,m,n] = (( 
			 l == 1  ? -(b[cnt,1] + b[cnt+1,1]) :
			 l == 2  ?   b[cnt,1] + b[cnt-1,1] - b[cnt+1,1] : 
			 l == Nx ?   b[cnt,1] : b[cnt,1]   - b[cnt+1,1])
			 +(
			 m == 1  ? -(b[cnt,2] + b[cnt+Nx,2]) :
			 m == 2  ?   b[cnt,2] + b[cnt-Nx,2] - b[cnt+Nx,2] : 
			 m == Ny ?   b[cnt,2] : b[cnt,2]    - b[cnt+Nx,2] )
			 +(
			 n == 1  ? -(b[cnt,3] + b[cnt+Nxy,3]) :
			 n == 2  ?   b[cnt,3] + b[cnt-Nxy,3] - b[cnt+Nxy,3] : 
			 n == Nz ?   b[cnt,3] : b[cnt,3]     - b[cnt+Nxy,3] ))
		end

	end
end
