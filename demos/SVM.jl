using RegLS

srand(123)
Na,Nb = 200,1500
a = [0+randn(Na)  3+ 4+4*randn(Na)]    
b = [2+randn(Nb)  3+-1+randn(Nb)]

X = OptVar(zeros(3))
A = [[a[:,1];b[:,1]] [a[:,2];b[:,2]] [ones(Na);ones(Nb)]]

bi = [ones(Na);-ones(Nb)]
d = 20*ones(3)  # regularization
d[end] .*= 1e-1 # reduce regularization in intercept

@time x,slv = minimize(ls(diagop(X,d))+hingeloss(A*X,bi), ZeroFPR(tol = 1e-8, verbose = 0))
show(slv)
@time x,slv = minimize(ls(diagop(X,d))+hingeloss(A*X,bi), FPG(tol = 1e-8, verbose = 0))
show(slv)
xx = linspace(-30,30,10)
x2 = copy(x)

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


