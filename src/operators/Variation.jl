export Variation

immutable Variation{N} <: LinearOperator
	domainType::Type
	dim_in::Tuple
	dim_out::Tuple
	Variation(domainType, dim_in) = new(domainType, dim_in, (prod(dim_in), length(dim_in)) )
end

size(L::Variation) = (L.dim_out, L.dim_in)

# Constructors

Variation(domainType::Type, dim_in::Tuple) = Variation{length(dim_in),dir}(domainType, dim_in)
Variation(domainType::Type, dim_in::Vararg{Int}) = Variation(domainType, dim_in)

Variation(dim_in::Tuple) = Variation{length(dim_in)}(Float64, dim_in)
Variation(dim_in::Vararg{Int}) = Variation(dim_in)

Variation{T,N}(x::AbstractArray{T,N})  = Variation{N}(eltype(x), size(x) ) 


# Operators

function A_mul_B!{T}(y::AbstractArray{T,2}, A::Variation{2}, b::AbstractArray{T,2})
	cnt = 0
	for m = 1:size(b,2), l = 1:size(b,1) 
		cnt += 1
		y[cnt,1] = l == 1 ? b[l+1,m]-b[l,m] : b[l,m]-b[l-1,m]
		y[cnt,2] = m == 1 ? b[l,m+1]-b[l,m] : b[l,m]-b[l,m-1]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,2}, A::Variation{2}, b::AbstractArray{T,2})
		
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

function A_mul_B!{T}(y::AbstractArray{T,2},A::Variation{3},b::AbstractArray{T,3})
	cnt = 0
	for n = 1:size(b,3), m = 1:size(b,2), l = 1:size(b,1) 
		cnt += 1
		y[cnt,1] = l == 1 ? b[l+1,m,n]-b[l,m,n] : b[l,m,n]-b[l-1,m,n]
		y[cnt,2] = m == 1 ? b[l,m+1,n]-b[l,m,n] : b[l,m,n]-b[l,m-1,n]
		y[cnt,3] = n == 1 ? b[l,m,n+1]-b[l,m,n] : b[l,m,n]-b[l,m,n-1]
	end
end

function Ac_mul_B!{T}(y::AbstractArray{T,3},A::Variation{3},b::AbstractArray{T,2})
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

# Properties

fun_name(L::Variation)  = "Variation Operator"
