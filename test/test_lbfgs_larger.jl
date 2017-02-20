@printf("Testing L-BFGS routines on a larger problem\n")

n = 10000
dens = 100/n

Q = sprandn(n, n, dens/2)
Q = 0.5*(Q+Q') + spdiagm(ones(n), 0) # this is spd
q = randn(n)

N = 15 # number of steps

mem = 5;

x = OptVar(zeros(n))
A = lbfgs(x,mem)
show(A)

x_old = 0;
grad_old = 0;
dir  = zeros(n) # matrix of directions (to be filled in)

for i = 1:N
    x = randn(n)
    grad = Q*x + q
    if i > 1
	    @time update!(A, x, x_old, grad, grad_old)
			A_mul_B!(dir,A,grad)
    else
	    copy!(dir, -grad)
    end
    @assert vecdot(dir, grad) < 0
    x_old = x
    grad_old = grad
end
