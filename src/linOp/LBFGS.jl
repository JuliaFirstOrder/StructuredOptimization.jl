export lbfgs, update!

immutable LBFGS{D1,D2} <:LinearOp{D1,D2}
	x::OptVar
  mem::Int64
	currmem::Array{Int64,1}
	curridx::Array{Int64,1}
  s::AbstractArray
  y::AbstractArray
  s_m::Array
  y_m::Array
  ys_m::Array{Float64,1}
  alphas::Array{Float64,1}
	H::Array{Float64,1}
end
size(A::LBFGS) = (size(A.x),size(A.x))
fun_name(A::LBFGS)  = "LBFGS Operator"

lbfgs(x::AbstractArray, args...) = lbfgs(OptVar(x), args...)

function lbfgs{D1}(x::OptVar{D1}, 
									 mem::Int64, 
									 currmem::Array{Int64,1},
									 curridx::Array{Int64,1},
									 ys_m::Array, 
									 alphas::Array, 
									 H::Array)
	s_m = Array{Array{D1},1}(mem)
  y_m = Array{Array{D1},1}(mem)
	for i in eachindex(s_m)
		s_m[i] = similar(x.x)
		y_m[i] = similar(x.x)
	end
	s = similar(x.x)
	y = similar(x.x)
	return LBFGS{D1,D1}(x, mem, currmem, curridx, s, y, s_m, y_m, ys_m, alphas, H)
end

function lbfgs{D1}(x::OptVar{D1}, mem::Int64)
	ys_m = zeros(Float64,mem)
	alphas = zeros(Float64,mem)
	lbfgs(x, mem, [0], [0], ys_m, alphas, [1.])
end

#create an array of LBFGS Op. which all share some stuff
function lbfgs{T<:AbstractArray}(x::Array{T,1},mem::Int64)
	LBFGS_col = Array{LBFGS,1}(length(x))
	LBFGS_col[1] = lbfgs(x[1],mem)
	for i = 2:length(x)
		LBFGS_col[i] = lbfgs(x[i],mem,
											   LBFGS_col[1].currmem,
												 LBFGS_col[1].curridx,
												 LBFGS_col[1].ys_m,
												 LBFGS_col[1].alphas,
												 LBFGS_col[1].H
												 ) 
	end
	return LBFGS_col
end

function update!{D1,D2}(A::LBFGS{D1,D2}, x::Array, x_prev::Array, gradx::Array, gradx_prev::Array)

	A.s .= (-).(x, x_prev)
	A.y .= (-).(gradx, gradx_prev)
	ys = real(vecdot(A.s,A.y))

	if ys > 0
		A.curridx[1] += 1
		if A.curridx[1] > A.mem A.curridx[1] = 1 end
		A.currmem[1] += 1
		if A.currmem[1] > A.mem A.currmem[1] = A.mem end

		copy!(A.s_m[A.curridx[1]], A.s)
		copy!(A.y_m[A.curridx[1]], A.y)
		A.ys_m[A.curridx[1]] = ys
		A.H[1] = ys/real(vecdot(A.y,A.y))
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

function A_mul_B!(d::AbstractArray, A::LBFGS, gradx::AbstractArray)
	d .= (-).(gradx)
	idx = A.curridx[1]
	for i=1:A.currmem[1]
		A.alphas[idx] = real(vecdot(A.s_m[idx], d))/A.ys_m[idx]
		d .= (-).(d, (*).(A.alphas[idx], A.y_m[idx]))
		idx -= 1
		if idx == 0 idx = A.mem end
	end
	d .= (*).(A.H[1], d)
	for i=1:A.currmem[1]
		idx += 1
		if idx > A.mem idx = 1 end
		beta = real(vecdot(A.y_m[idx], d))/A.ys_m[idx]
		d .= (+).(d, (*).((A.alphas[idx]-beta), A.s_m[idx]))
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
