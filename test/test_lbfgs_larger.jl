@printf("Testing L-BFGS routines on a larger problem\n")

n = 10000
dens = 100/n

Q = sprandn(n, n, dens/2)
Q = 0.5*(Q+Q') + spdiagm(ones(n), 0) # this is spd
q = randn(n)

N = 15 # number of steps

mem = 5;

lbfgs = RegLS.LBFGS(mem,zeros(n))
x_old = 0;
grad_old = 0;

for i = 1:N
    x = randn(n)
    grad = Q*x + q
    if i > 1
	    @time RegLS.push!(lbfgs, x, x_old, grad, grad_old)
    else
	    @time RegLS.push!(lbfgs, grad)
    end
    dir = lbfgs.d
    @assert vecdot(dir, grad) < 0
    x_old = x
    grad_old = grad
end
