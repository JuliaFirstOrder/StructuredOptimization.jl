module LBFGS

type Storage
  mem::Int64
  currmem::Int64
  curridx::Int64
  s_m::Array{Array,1}
  y_m::Array{Array,1}
  ys_m::Array{Float64,1}
end

function create(mem::Int64)
	s_m = Array{Array,1}(mem)
	y_m = Array{Array,1}(mem)
	ys_m = zeros(mem)
	Storage(mem, 0, 0, s_m, y_m, ys_m)
end

function push(obj::Storage, s::Array{Float64}, y::Array{Float64})
	obj.curridx += 1
  if obj.curridx > obj.mem obj.curridx = 1 end
	obj.currmem += 1
  if obj.currmem > obj.mem obj.currmem = obj.mem end
	obj.s_m[obj.curridx] = s
	obj.y_m[obj.curridx] = y
	obj.ys_m[obj.curridx] = vecdot(s,y)
end

function matvec(obj::Storage, H::Float64, g::Array{Float64})
	# TODO maybe this operation alone can be taken from libLBFGS
	q = g
	alphas = zeros(obj.mem)
	idx = obj.curridx
	for i=1:obj.currmem
		alphas[idx] = vecdot(obj.s_m[idx],q)/obj.ys_m[idx]
		q = q - alphas[idx]*obj.y_m[idx]
		idx = idx - 1
		if idx == 0 idx = obj.mem end
	end
	z = H*q # modify here
	for i=1:obj.currmem
		idx = idx + 1
		if idx > obj.mem idx = 1 end
		beta = vecdot(obj.y_m[idx],z)/obj.ys_m[idx]
		z = z + (alphas[idx]-beta)*obj.s_m[idx]
	end
	return z
end

function reset(obj::Storage)
	obj.currmem = 0
	obj.curridx = 0
end

end
