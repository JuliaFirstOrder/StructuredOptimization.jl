using ProximalOperators
using RegLS

srand(123)
Na,Nb = 200,1500
a = [0+randn(Na)  20+randn(Na)]    #working
b = [-20+randn(Nb)  -10+randn(Nb)]
a = [0+randn(Na)  20+randn(Na)]    # not working
b = [0+randn(Nb)  -10+randn(Nb)]

#Na,Nb = 1,1
#a = [0 10]    #working
#b = [0 5]
#a = [0 10]    #not working
#b = [0 -30]
#
x = zeros(3)
A = [[a[:,1];b[:,1]] [a[:,2];b[:,2]] [ones(Na);ones(Nb)]]

g = HingeLoss([ones(Na);-ones(Nb)],1.)

@time x,slv = solve(zeros(x), g, A, zeros(Na+Nb), ZeroFPR(tol = 1e-8))
show(slv)
@time x,slv = solve(zeros(x), g, A, zeros(Na+Nb), FPG(tol = 1e-8))
show(slv)
xx = linspace(-30,30,10)
x2 = copy(x)

L = 1/norm(A)^2
tau = 1/norm(A) 
sigma = 1/norm(A) 
theta = 1. 
x = zeros(3)
xbar = zeros(3)
z = zeros(Na+Nb)
for k = 1:10000

	z += sigma*(A*xbar) 
	prox!(Conjugate(g),z,sigma)
	xprev = copy(x)
	prox!(SqrNormL2([1.;1.;1e-4]),x-tau*A'*z,x) 
	xbar = x+theta*(x-xprev)

end


theta = atan(x[2]/x[1])
marginy = cos(theta)/norm(x[1:2])
marginx = sin(theta)/norm(x[1:2])

using PyPlot
figure()
plot(a[:,1],a[:,2], "r*")
plot(b[:,1],b[:,2], "ko")
plot(xx,-(x[1]*xx+x[3])./x[2])
plot(xx+marginy,-(x[1]*xx+x[3])./x[2]+marginx)
plot(xx-marginy,-(x[1]*xx+x[3])./x[2]-marginx )
xlim([-30;30])
ylim([-30;30])


