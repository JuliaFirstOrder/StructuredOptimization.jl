export LBFGS, update!

immutable LBFGS{N,M, T<:RealOrComplex, A<:AbstractArray{T}} <:LinearOperator
	domainType::Type
	dim::NTuple{N,Int}
	currmem::Array{Int64,1}
	curridx::Array{Int64,1}
	s::A
	y::A
	s_m::NTuple{M,A}
	y_m::NTuple{M,A}
	ys_m::Array{T,1}
	alphas::Array{T,1}
	H::Array{T,1}
end
size(A::LBFGS) = (A.dim,A.dim)

function LBFGS{T<:RealOrComplex, A<:AbstractArray{T}}(x::A, mem::Int)

	s_m = ([similar(x) for i = 1:mem]...)
	y_m = ([similar(x) for i = 1:mem]...)

	s = similar(x)
	y = similar(x)

	ys_m = zeros(Float64,mem)
	alphas = zeros(Float64,mem)

	return LBFGS{ndims(x),mem,T,A}(eltype(x), size(x), [0], [0], s, y, s_m, y_m, ys_m, alphas, [1.])
end

#create an array of LBFGS Op. which all share some stuff
function LBFGS{T<:AbstractArray}(x::Array{T,1},mem::Int64)
	LBFGS_col = Array{LBFGS,1}(length(x))
	LBFGS_col[1] = LBFGS(x[1],mem)
	for i = 2:length(x)
		LBFGS_col[i] = LBFGS(x[i],mem,
		                     LBFGS_col[1].currmem,
				     LBFGS_col[1].curridx,
				     LBFGS_col[1].ys_m,
				     LBFGS_col[1].alphas,
				     LBFGS_col[1].H) 
	end
	return LBFGS_col
end

function update!{N,M, 
		 T<:RealOrComplex, 
		 A<:AbstractArray{T}}(
			L::LBFGS{N,M,T,A}, 
			x::A, 
			x_prev::A, 
			gradx::A, 
			gradx_prev::A)

	L.s .= (-).(x, x_prev)
	L.y .= (-).(gradx, gradx_prev)
	ys = real(vecdot(L.s,L.y))

	if ys > 0
		L.curridx[1] += 1
		if L.curridx[1] > M L.curridx[1] = 1 end
		L.currmem[1] += 1
		if L.currmem[1] > M L.currmem[1] = M end

		L.s_m[L.curridx[1]] .=  L.s
		L.y_m[L.curridx[1]] .=  L.y
		L.ys_m[L.curridx[1]] = ys
		L.H[1] = ys/real(vecdot(L.y,L.y))
	end

end

function update!{T<:AbstractArray}(A::Array{LBFGS,1}, 
				   x::Array{T,1}, 
				   x_prev::Array{T,1}, 
				   gradx::Array{T,1}, 
				   gradx_prev::Array{T,1})

	ys = 0.
	for i in eachindex(A)
		A[i].s .= (-).(x[i], x_prev[i])
		A[i].y .= (-).(gradx[i], gradx_prev[i])
		ys += real(vecdot(A[i].s, A[i].y))
	end

	if ys > 0
		A[1].curridx[1] += 1
		if A[1].curridx[1] > A[1].mem A[1].curridx[1] = 1 end
		A[1].currmem[1] += 1
		if A[1].currmem[1] > A[1].mem A[1].currmem[1] = A[1].mem end

		A[1].H[1] = 0.
		for i in eachindex(A)
			copy!(A[i].s_m[A[1].curridx[1]], A[i].s)
			copy!(A[i].y_m[A[1].curridx[1]], A[i].y)
			A[1].H[1] += real(vecdot(A[i].y,A[i].y))
		end
		A[1].ys_m[A[1].curridx[1]] = ys
		A[1].H[1] = (A[1].H[1]/ys)^(-1)
	end

end

function A_mul_B!{N, M, T<:RealOrComplex, A<:AbstractArray{T}}(d::A, L::LBFGS{N,M,T,A}, gradx::A)
	d .= (-).(gradx)
	idx = L.curridx[1]
	for i=1:L.currmem[1]
		L.alphas[idx] = real(vecdot(L.s_m[idx], d))/L.ys_m[idx]
		d .= (-).(d, (*).(L.alphas[idx], L.y_m[idx]))
		idx -= 1
		if idx == 0 idx = M end
	end
	d .= (*).(L.H[1], d)
	for i=1:L.currmem[1]
		idx += 1
		if idx > M idx = 1 end
		beta = real(vecdot(L.y_m[idx], d))/L.ys_m[idx]
		d .= (+).(d, (*).((L.alphas[idx]-beta), L.s_m[idx]))
	end
end

function A_mul_B!{T<:AbstractArray}(d::Array{T,1}, A::Array{LBFGS,1}, gradx::Array{T,1})
	d .= (-).(gradx)
	idx = A[1].curridx[1]
	for i=1:A[1].currmem[1]
		
		A[1].alphas[idx] = 0.
		for i in eachindex(A)
			A[1].alphas[idx] += real(vecdot(A[i].s_m[idx], d[i]))
		end
		A[1].alphas[idx] /= A[1].ys_m[idx]
		
		for i in eachindex(A)
			d[i] .= (-).(d[i], (*).(A[1].alphas[idx], A[i].y_m[idx]))
		end
		
		idx -= 1
		if idx == 0 idx = A[1].mem end
	end
	d .= (*).(A[1].H[1], d)
	for i=1:A[1].currmem[1]
		idx += 1
		if idx > A[1].mem idx = 1 end
		
		beta = 0.
		for i in eachindex(A)
			beta += real(vecdot(A[i].y_m[idx], d[i]))
		end
		beta /= A[1].ys_m[idx]

		for i in eachindex(A)
			d[i] .= (+).(d[i], (*).((A[1].alphas[idx]-beta), A[i].s_m[idx]))
		end
	end
end

function reset(A::LBFGS)
	A.currmem[1] = 0
	A.curridx[1] = 0
end

function reset(A::Array{LBFGS,1})
	A[1].currmem[1] = 0
	A[1].curridx[1] = 0
end

fun_name(A::LBFGS)  = "LBFGS Operator"
