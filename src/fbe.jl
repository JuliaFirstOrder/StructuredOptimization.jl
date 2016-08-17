# This function is for debugging purposes

function fbe(x::Array{Float64}, gamma::Float64, L, Ladj, b, proxg)
	resx = L(x) - b
	fx = 0.5*vecnorm(resx)^2
	gradx = Ladj(resx)
	xbar, gxbar = proxg(x-gamma*gradx, gamma)
	r = x - xbar
	normr = vecnorm(r)
	return fx + gxbar - vecdot(gradx,r) + 1/(2*gamma)*normr^2
end
