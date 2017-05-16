@printf("Testing L-BFGS routines on a larger problem\n")

n = 10000
dens = 100/n

Q = sprandn(n, n, dens/2)
Q = 0.5*(Q+Q') + spdiagm(ones(n), 0) # this is spd
q = randn(n)

N = 15 # number of steps

mem = 5;

x = zeros(n)
A = LBFGS(x,mem)
println(A)
println()
nh = round(Int64,n/2)
x2 = (zeros(nh),zeros(nh))

A2 = LBFGS(x2,mem)
#println(A2)

x_old = randn(n)
x_old2 = (randn(nh),randn(nh))
grad_old = randn(n)
grad_old2 = (randn(nh),randn(nh))
grad = randn(n)
grad2 = (randn(nh),randn(nh))

dir  = zeros(n) # matrix of directions (to be filled in)
dir2 = (zeros(round(Int64,n/2)),zeros(round(Int64,n/2))) # matrix of directions (to be filled in)

for i = 1:N
    x = randn(n)
    x2 = (x[1:nh],x[nh+1:end])
    grad = Q*x + q
    grad2 = (grad[1:nh],grad[nh+1:end])
    if i > 1
	    update!(A,   x, x_old,   grad,  grad_old)
	    @time update!(A2, x2, x_old2, grad2, grad_old2)
	    A_mul_B!(dir,A,grad)
	    A_mul_B!(dir2,A2,grad2)
    else
	    copy!(dir, -grad)
	    dir2 = deepcopy((-).(grad2))
    end
    @assert vecdot(dir, grad) < 0
    x_old  = x
    x_old2 = deepcopy(x2)

    grad_old  = grad
    grad_old2 = deepcopy(grad2)
end

@test norm(dir[1:nh]-dir2[1])<1e-8
@test norm(dir[nh+1:end]-dir2[2])<1e-8

x2 = (randn(nh),randn(nh)+randn(nh)*im)
x_old2 = (randn(nh),randn(nh)+randn(nh)*im)
grad2 = (randn(nh),randn(nh)+randn(nh)*im)
grad_old2 = (randn(nh),randn(nh)+randn(nh)*im)
A2 = LBFGS(x2,mem)
#println(A2)
	    
@time update!(A2, x2, x_old2, grad2, grad_old2)
@time update!(A2, x2, x_old2, grad2, grad_old2)
@time update!(A2, x2, x_old2, grad2, grad_old2)

