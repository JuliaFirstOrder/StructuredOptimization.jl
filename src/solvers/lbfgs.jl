import Base: push!

type LBFGS{T <: AbstractArray}
  mem::Int64
  currmem::Int64
  curridx::Int64
  s::T
  y::T
  s_m::Array{T}
  y_m::Array{T}
  ys_m::Array{Float64,1}
  alphas::Array{Float64,1}
  H::Float64
  d::T
end

function LBFGS{T <: AbstractArray}(mem::Int64,x::T)
	ys_m = zeros(Float64,mem)
	alphas = zeros(Float64,mem)
  s_m = Array(T, mem)
  y_m = Array(T, mem)
	for i in eachindex(s_m)
		s_m[i] = deepsimilar(x)
		y_m[i] = deepsimilar(x)
	end
	s = deepsimilar(x)
	y = deepsimilar(x)
  d = deepcopy(x)
	LBFGS{T}(mem, 0, 0, s, y, s_m, y_m, ys_m, alphas, 1., d)
end

function push!(obj::LBFGS, gradx::Array)
	obj.d .= (-).(gradx)
end

function push!(obj::LBFGS, x::Array, x_prev::Array, gradx::Array, gradx_prev::Array)

	obj.s .= (-).(x, x_prev)
	obj.y .= (-).(gradx, gradx_prev)
	ys = real(deepvecdot(obj.s,obj.y))

	if ys > 0
		obj.curridx += 1
		if obj.curridx > obj.mem obj.curridx = 1 end
		obj.currmem += 1
		if obj.currmem > obj.mem obj.currmem = obj.mem end

		deepcopy!(obj.s_m[obj.curridx], obj.s)
		deepcopy!(obj.y_m[obj.curridx], obj.y)
		obj.ys_m[obj.curridx] = ys
		obj.H = ys/real(deepvecdot(obj.y,obj.y))
	end

	twoloop(obj, gradx)
end

function twoloop(obj::LBFGS, gradx::Array)
	obj.d .= (-).(gradx)
	idx = obj.curridx
	for i=1:obj.currmem
		obj.alphas[idx] = real(deepvecdot(obj.s_m[idx], obj.d))/obj.ys_m[idx]
		obj.d .= (-).(obj.d, (*).(obj.alphas[idx], obj.y_m[idx]))
		idx -= 1
		if idx == 0 idx = obj.mem end
	end
	obj.d .= (*).(obj.H, obj.d)
	for i=1:obj.currmem
		idx += 1
		if idx > obj.mem idx = 1 end
		beta = real(deepvecdot(obj.y_m[idx], obj.d))/obj.ys_m[idx]
		obj.d .= (+).(obj.d, (*).((obj.alphas[idx]-beta), obj.s_m[idx]))
	end
end

function reset(obj::LBFGS)
	obj.currmem = 0
	obj.curridx = 0
end
