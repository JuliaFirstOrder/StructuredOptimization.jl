module LBFGS

type Storage
  mem::Int64
  currmem::Int64
  curridx::Int64
  s_m::Array{Array,1}
  y_m::Array{Array,1}
  ys_m::Array{Float64,1}
  alphas::Array{Float64,1}
  H::Float64
  d::Array
end

function create(mem::Int64,x::Array)
	s_m = Array{Array,1}(mem)
	y_m = Array{Array,1}(mem)
	ys_m = zeros(Float64,mem)
	alphas = zeros(Float64,mem)
	for i = 1:mem  s_m[i] = zeros(x); y_m[i] = zeros(x); end
	Storage(mem, 0, 0, s_m, y_m, ys_m, alphas, 1., zeros(x))
end

function push!(obj::Storage, gradx::Array)
	obj.d[:] = gradx
end

function push!(obj::Storage, x::Array, x_prev::Array, gradx::Array, gradx_prev::Array)
	
	obj.curridx += 1
	if obj.curridx > obj.mem obj.curridx = 1 end
	obj.currmem += 1
	if obj.currmem > obj.mem obj.currmem = obj.mem end

	obj.s_m[obj.curridx][:] = x-x_prev
	obj.y_m[obj.curridx][:] = gradx-gradx_prev
	obj.ys_m[obj.curridx] = real(vecdot(obj.s_m[obj.curridx], obj.y_m[obj.curridx]))
	obj.H = obj.ys_m[obj.curridx]/real(vecdot(obj.y_m[obj.curridx],obj.y_m[obj.curridx]))
	obj.d[:] = matvec(obj,gradx)
end

function matvec(obj::Storage, gradx::Array)
	# TODO maybe this operation alone can be taken from libLBFGS
	d = -gradx
	idx = obj.curridx
	for i=1:obj.currmem
		obj.alphas[idx] = real(vecdot(obj.s_m[idx],d))/obj.ys_m[idx]
		d -= obj.alphas[idx]*obj.y_m[idx]
		idx -= 1
		if idx == 0 idx = obj.mem end
	end
	d = obj.H*d # modify here
	for i=1:obj.currmem
		idx += 1
		if idx > obj.mem idx = 1 end
		beta = real(vecdot(obj.y_m[idx],d))/obj.ys_m[idx]
		d += (obj.alphas[idx]-beta)*obj.s_m[idx]
	end
	return d
end

function reset(obj::Storage)
	obj.currmem = 0
	obj.curridx = 0
end

end
