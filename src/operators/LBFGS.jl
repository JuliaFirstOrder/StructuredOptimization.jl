export LBFGS, update!

#TODO make Ac_mul_B!

type LBFGS{M, T<:RealOrComplex, A <: AbstractArray{T}, N} <:LinearOperator
	domainType::Type
	dim::NTuple{N,Int}
	currmem::Int
	curridx::Int
	s::A
	y::A
	s_m::NTuple{M,A}
	y_m::NTuple{M,A}
	ys_m::Array{Float64,1}
	alphas::Array{Float64,1}
	H::Float64
end

# Constructors

function LBFGS{T<:RealOrComplex, A<:AbstractArray{T}}(x::A, mem::Int)

	s_m = ([similar(x) for i = 1:mem]...)
	y_m = ([similar(x) for i = 1:mem]...)

	s = similar(x)
	y = similar(x)

	ys_m = zeros(Float64,mem)
	alphas = zeros(Float64,mem)

	return LBFGS{mem,T,A,ndims(x)}(eltype(x), size(x), 0, 0, s, y, s_m, y_m, ys_m, alphas, 1.)
end

function update!{M,
		 T<:RealOrComplex,
		 A<:AbstractArray{T},N}(
			L::LBFGS{M,T,A,N},
			x::A,
			x_prev::A,
			gradx::A,
			gradx_prev::A)


	ys = update_s_y(L,x,x_prev,gradx,gradx_prev)

	if ys > 0
		L.curridx += 1
		if L.curridx > M L.curridx = 1 end
		L.currmem += 1
		if L.currmem > M L.currmem = M end


		yty = update_s_m_y_m(L,L.curridx)
		L.ys_m[L.curridx] = ys
		L.H = ys/yty
	end
	return L
end

function update_s_y{M,T,A,N}(L::LBFGS{M,T,A,N}, x::A, x_prev::A, gradx::A, gradx_prev::A)
	L.s .= (-).(x, x_prev)
	L.y .= (-).(gradx, gradx_prev)
	ys = real(vecdot(L.s,L.y))
	return ys
end

function update_s_m_y_m{M,T,A,N}(L::LBFGS{M,T,A,N}, curridx::Int)
	L.s_m[curridx] .=  L.s
	L.y_m[curridx] .=  L.y

	yty = real(vecdot(L.y,L.y))
	return yty
end

function A_mul_B!{M, T<:RealOrComplex, A<:AbstractArray{T},N}(d::A, L::LBFGS{M,T,A,N}, gradx::A)
	d .= (-).(gradx)
	idx = loop1!(d,L)
	d .= (*).(L.H, d)
	d = loop2!(d,idx,L)
end

function loop1!{M,T,A,N}(d::A, L::LBFGS{M,T,A,N})
	idx = L.curridx
	for i=1:L.currmem
		L.alphas[idx] = real(vecdot(L.s_m[idx], d))/L.ys_m[idx]
		d .-= L.alphas[idx].*L.y_m[idx]
		idx -= 1
		if idx == 0 idx = M end
	end
	return idx
end

function loop2!{M,T,A,N}(d::A, idx::Int, L::LBFGS{M,T,A,N})
	for i=1:L.currmem
		idx += 1
		if idx > M idx = 1 end
		beta = real(vecdot(L.y_m[idx], d))/L.ys_m[idx]
		d .+=  (L.alphas[idx].-beta).*L.s_m[idx]
	end
	return d
end

# Properties

size(A::LBFGS) = (A.dim,A.dim)

fun_name(A::LBFGS)  = "LBFGS Operator"
