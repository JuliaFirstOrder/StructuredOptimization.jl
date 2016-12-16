using RegLS
using ProximalOperators
using Base.Test

n = 2000 # number of features
m = 300 # number of data points

w = sprandn(n, 0.3) # N(0,1), 30% sparse
v = randn() # random intercept

X = sprandn(m, n, 10/n)
btrue = sign(X*w + v)

b = sign(X*w + v + sqrt(0.1)*randn(m)); # labels with noise

A = [X ones(m)]; # add bias term coefficient

ratio = sum(b .== 1)/m
idx_hi = (b .== 1)
idx_lo = (b .== -1)
lam = 0.1 * norm((1-ratio)*sum(A[idx_hi, :], 1) + ratio*sum(A[idx_lo, :], 1), Inf)

y, slv = solve(zeros(n+1), HingeLoss(b, 1/lam), A)
