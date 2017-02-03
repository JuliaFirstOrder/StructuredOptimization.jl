using ProximalOperators
using RegLS

srand(123)
Na,Nb = 200,1500
a = [0+randn(Na)  20+randn(Na)]    #working
b = [-20+randn(Nb)  -10+randn(Nb)]
#a = [0+randn(Na)  20+randn(Na)]    # not working
#b = [0+randn(Nb)  -10+randn(Nb)]

#Na,Nb = 1,1
#a = [0 10]    #working
#b = [0 5]
#a = [0 10]    #not working
#b = [0 -15]

x = zeros(3)
A = [[a[:,1];b[:,1]] [a[:,2];b[:,2]] [ones(Na);ones(Nb)]]

g = HingeLoss([ones(Na);-ones(Nb)],1e9)

@time x, = solve(zeros(x), g, A, zeros(Na+Nb), ZeroFPR(tol = 1e-8))
xx = linspace(-30,30,10)

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
#

