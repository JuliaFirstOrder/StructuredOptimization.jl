export lbfgs, update!

type LBFGS{D1,D2} <:LinearOp{D1,D2}
	x::OptVar
  mem::Int64
  currmem::Int64
  curridx::Int64
  s::AbstractArray
  y::AbstractArray
  s_m::Array
  y_m::Array
  ys_m::Array{Float64,1}
  alphas::Array{Float64,1}
  H::Float64
end
size(A::LBFGS) = (size(A.x),size(A.x))
fun_name(A::LBFGS)  = "LBFGS Operator"

function lbfgs{D1}(x::OptVar{D1}, mem::Int64)
	ys_m = zeros(Float64,mem)
	alphas = zeros(Float64,mem)
	s_m = Array{Array{D1},1}(mem)
  y_m = Array{Array{D1},1}(mem)
	for i in eachindex(s_m)
		s_m[i] = similar(x.x)
		y_m[i] = similar(x.x)
	end
	s = similar(x.x)
	y = similar(x.x)
	return LBFGS{D1,D1}(x, mem, 0, 0, s, y, s_m, y_m, ys_m, alphas, 1.)
end

function lbfgs(xyz::Array{OptVar,1},mem::Int64)
	LBFGSs = Array{LBFGS,1}(length(xyz))
	for i in eachindex(LBFGSs)
		LBFGSs[i] = lbfgs(xyz[i],mem)
	end
	return LBFGSs
end

function update!{D1,D2}(A::LBFGS{D1,D2}, x::Array, x_prev::Array, gradx::Array, gradx_prev::Array)

	A.s .= (-).(x, x_prev)
	A.y .= (-).(gradx, gradx_prev)
	ys = real(deepvecdot(A.s,A.y))

	if ys > 0
		A.curridx += 1
		if A.curridx > A.mem A.curridx = 1 end
		A.currmem += 1
		if A.currmem > A.mem A.currmem = A.mem end

		deepcopy!(A.s_m[A.curridx], A.s)
		deepcopy!(A.y_m[A.curridx], A.y)
		A.ys_m[A.curridx] = ys
		A.H = ys/real(deepvecdot(A.y,A.y))
	end

end

function update!{T<:AbstractArray}(As::Array{LBFGS,1}, 
																	 x::Array{T,1}, 
																	 x_prev::Array{T,1}, 
																	 gradx::Array{T,1}, 
																	 gradx_prev::Array{T,1})

	for i = 1:length(As)
		update!(As[i],x[i],x_prev[i],gradx[i],gradx_prev[i])
	end
end

function A_mul_B!(d::AbstractArray, A::LBFGS, gradx::AbstractArray)
	d .= (-).(gradx)
	idx = A.curridx
	for i=1:A.currmem
		A.alphas[idx] = real(deepvecdot(A.s_m[idx], d))/A.ys_m[idx]
		d .= (-).(d, (*).(A.alphas[idx], A.y_m[idx]))
		idx -= 1
		if idx == 0 idx = A.mem end
	end
	d .= (*).(A.H, d)
	for i=1:A.currmem
		idx += 1
		if idx > A.mem idx = 1 end
		beta = real(deepvecdot(A.y_m[idx], d))/A.ys_m[idx]
		d .= (+).(d, (*).((A.alphas[idx]-beta), A.s_m[idx]))
	end
end

function A_mul_B!{T<:AbstractArray}(d::Array{T,1},As::Array{LBFGS,1},gradx::Array{T,1})
	for i = 1:length(As)
		A_mul_B!(d[i],As[i],gradx[i])
	end
end


function reset(A::LBFGS)
	A.currmem = 0
	A.curridx = 0
end
